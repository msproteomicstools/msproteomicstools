#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013, ETH Zurich
For a full list of authors, refer to the file AUTHORS.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""

from sys import stdout
import numpy
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
import msproteomicstoolslib.data_structures.PeakGroup
import msproteomicstoolslib.math.Smoothing as smoothing

def getDistanceMatrix(exp, multipeptides, initial_alignment_cutoff):
    """Compute distance matrix of all runs.

    Computes a n x n distance matrix between all runs of an experiment. The
    reported distance is the standard deviation of the aligned values when
    fitted using a linear regression fit.

    Args:
        exp(MRExperiment): a collection of runs
        multipeptides(list(Multipeptide)): a list of
            multipeptides containing the matching of precursors across runs.
        initial_alignment_cutoff(float): a filtering cutoff (in q-value) to
            specify which points should be used for the calculation of the
            distance. In general, only identification which are very certain
            should be used for this and a q-value of 0.0001 is recommended --
            given that there are enough points.

    Returns:
        None
        numpy (n x n) matrix(float): distance matrix
    """

    spl_aligner = SplineAligner(initial_alignment_cutoff)
    dist_matrix = numpy.zeros(shape=(len(exp.runs),len(exp.runs)))
    for i in range(len(exp.runs)):
        for j in range(len(exp.runs)):
            if i == j:
                dist_matrix[i,j] = 0
                continue 

            idata, jdata = spl_aligner._getRTData(exp.runs[i], exp.runs[j], multipeptides)

            if len(idata) == 0:
                dist_matrix[i,j] = 0
                continue

            # Get linear alignment
            smlin = smoothing.SmoothingLinear()
            smlin.initialize(idata, jdata)
            idata_lin_aligned = smlin.predict(idata)

            # Use stdev to estimate distance between two runs
            stdev_lin = numpy.std(numpy.array(jdata) - numpy.array(idata_lin_aligned))
            dist_matrix[i,j] = stdev_lin

    return dist_matrix

class TreeConsensusAlignment():
    """ Multiple run alignment using a minimum spanning tree (MST).

    This class will align features across multiple runs using a strictly local
    approach. It uses a minimum spanning tree (MST) as input which is expected
    to allow traversal of all runs, only connecting the most similar runs. This
    should allow accurate local alignment between two very similar runs and
    transfer of identification with high accuracy. Specifically, this should
    for scalability to a large number of dissimilar runs where approaches that
    rely on a single reference run for alignment might give less accurate
    results.

    Briefly, the algorithm will choose the best scoring peakgroup as a seed and
    start to traverse the MST from this seed. At each node, it will add the
    best matching peakgroup (by score, within a specified retention time
    window) to the result. After traversing all nodes, a new seed can be chosen
    among the peakgroups not yet belonging to a cluster and the process can be
    repeated to produce multiple clusters.

    For example, consider a case of 5 LC-MS/MS runs where 6 different feature
    (peakgroups) were found in each run (not all peakgroups were found in all
    runs):

    Run 1:
    ---------- pg1_1 ------- pg1_2 --- pg1_3 - pg1_4 --------- pg1_5 ------- pg1_6

    Run 2:
    ----------- pg2_1 ------- pg2_2 --- pg2_3 - pg2_4 --------- pg2_5 ------ pg2_6

    Run 3:
    ----- pg3_0 ---------- pg3_1 ---- pg3_2 - pg3_3 - pg3_4 ----- pg3_5 ---- pg3_6

    Run 4:
    ----------- pg4_0 ---------- pg4_1 ------- pg4_2 --- pg4_3 - pg4_4 --------- pg4_5

    Run 5:
    -- pg5_0 ---------- pg5_1 ------- pg5_2 --- pg5_3 - pg5_4 --------- pg5_5 


    Assume that the corresponding MST looks like this:

                               /-- Run4
        Run1 -- Run2 -- Run3 --
                               \-- Run5

    This is a case where Run1 and Run2 are very similar and Run3 and Run4 are
    rather similar and should be easy to align. The algorithm will start with
    the "best" peakgroup overall (having the best probability score), assume
    this peakgroups is pg1_1 from Run 1. The algorithm will then use the
    alignment Run1-Run2 to infer that pg2_1 is the same signal as pg1_1 and add
    it to the group. Specifically, it will select the highest-scoring peakgroup
    within a narrow RT-window (max_rt_diff) in Run2 - note that if the
    RT-window is too wide, there is a certain chance of mis-matching, e.g. pg_2
    will be selected instead of pg2_1.  The alignment Run2-Run3 will be used to
    add pg3_1. Then a bifurcation in the tree occurs and Run3-Run4 as well as
    Run3-Run5 will be used to infer the identity of pg4_1 and pg5_1 and add
    them to the cluster.  In the end, the algorithm will report (pg1_1, pg2_1,
    pg3_1, pg4_1, pg5_1) as a consistent cluster across multiple runs. This
    process can be repeated with the next best peakgroup that is not yet part
    of a cluster (e.g. pg1_2) until no more peakgroups are left (no more
    peakgroups having a score below fdr_cutoff).

    Note how the algorithm only used binary alignments and purely local
    alignments of the runs that are most close to each other. This stands in
    contrast to approaches where a single reference is picked and then used for
    alignment which might align runs that are substantially different. On the
    other hand, a single error at one edge in the tree will propagate itself
    and could lead to whole subtrees that are wrongly aligned.
    """

    def __init__(self, max_rt_diff, fdr_cutoff, aligned_fdr_cutoff,
                 correctRT_using_pg=False, stdev_max_rt_per_run=None, verbose=False):
        """ Initialization with parameters

        Args:
            max_rt_diff(float): maximal difference in retention time to be used
                to look for a matching peakgroup in an adjacent run
            fdr_cutoff(float): maximal FDR that at least one peakgroup needs to
                reach (seed FDR)
            aligned_fdr_cutoff(float): maximal FDR that a peakgroup needs to
                reach to be considered for extension (extension FDR)
            correctRT_using_pg(bool): use the apex of the aligned peak group
                as the input for the next alignment during MST traversal
                (opposed to using the transformed RT plain)
            stdev_max_rt_per_run(float): use a different maximal RT tolerance
                for each alignment, depending on the goodness of the alignment.
                The RT tolerance used by the algorithm will be the standard
                deviation times stdev_max_rt_per_run.
        """

        self._max_rt_diff = max_rt_diff
        self._aligned_fdr_cutoff = aligned_fdr_cutoff
        self._fdr_cutoff = fdr_cutoff
        self._correctRT_using_pg = correctRT_using_pg
        self._stdev_max_rt_per_run = stdev_max_rt_per_run
        self.verbose = verbose

    def alignBestCluster(self, multipeptides, tree, tr_data):
        """Use the MST to report the first cluster containing the best peptide (overall).

        The algorithm will go through all multipeptides and mark those
        peakgroups which it deems to belong to the best peakgroup cluster (only
        the first cluster will be reported).

        Args:
            multipeptides(list(Multipeptide)): a list of
                multipeptides on which the alignment should be performed. After
                alignment, each peakgroup that should be quantified can be
                retrieved by calling get_selected_peakgroups() on the multipeptide.
            tree(list(tuple)): a minimum spanning tree (MST) represented as
                list of edges (for example [('0', '1'), ('1', '2')] ). Node names
                need to correspond to run ids.
            tr_data(format.TransformationCollection.LightTransformationData):
                structure to hold binary transformations between two different
                retention time spaces

        Returns:
            None
        """

        for m in multipeptides:

            # Find the overall best peptide
            best = m.find_best_peptide_pg()
            if self.verbose: 
                print "00000000000000000000000000000000000 new peptide (cluster)", m.get_peptides()[0].get_id()
                print " Best", best.print_out(), "from run", best.peptide.run.get_id()

            # Use this peptide to generate a cluster
            for pg_ in self._findAllPGForSeed(tree, tr_data, m, best, {}):
                pg_.select_this_peakgroup()

    def alignAllCluster(self, multipeptides, tree, tr_data):
        """Use the MST to report all cluster.

        Briefly, the algorithm will choose the best scoring peakgroup as a seed and
        start to traverse the MST from this seed. At each node, it will add the
        best matching peakgroup (by score, within a specified retention time
        window) to the result. After traversing all nodes, a new seed can be chosen
        among the peakgroups not yet belonging to a cluster and the process can be
        repeated to produce multiple clusters.  It will add clusters until no
        more peptides with an fdr score better than self._fdr_cutoff are left.

        Args:
            multipeptides(list(Multipeptide)): a list of
                multipeptides on which the alignment should be performed. After
                alignment, each peakgroup that should be quantified can be
                retrieved by calling get_selected_peakgroups() on the multipeptide.
            tree(list(tuple)): a minimum spanning tree (MST) represented as
                list of edges (for example [('0', '1'), ('1', '2')] ). Node names
                need to correspond to run ids.
            tr_data(format.TransformationCollection.LightTransformationData):
                structure to hold binary transformations between two different
                retention time spaces

        Returns:
            None
        """

        from msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm import Cluster
        for m in multipeptides:
            if self.verbose: 
                print "00000000000000000000000000000000000 new peptide (cluster)", m.get_peptides()[0].get_id()

            last_cluster = []
            already_seen = set([])
            stillLeft = [b for a in m.get_peptides() for b in a.get_all_peakgroups() 
                         if b.get_feature_id() + b.peptide.get_id() not in already_seen 
                         and b.get_fdr_score() < self._fdr_cutoff]
            clusters = []
            while len(stillLeft) > 0:
                best = min(stillLeft, key=lambda x: float(x.get_fdr_score()))
                if self.verbose: 
                    print "111111111111111111111111 new cluster"
                    print " Best", best.print_out(), "from run", best.peptide.run.get_id()

                last_cluster = self._findAllPGForSeed(tree, tr_data, m, best, already_seen)
                already_seen.update( set([ b.get_feature_id() + b.peptide.get_id() 
                                          for b in last_cluster if b is not None]) )
                stillLeft = [b for a in m.get_peptides() for b in a.get_all_peakgroups() 
                             if b.get_feature_id() + b.peptide.get_id() not in already_seen 
                             and b.get_fdr_score() < self._fdr_cutoff]
                clusters.append(Cluster(last_cluster))

            # select the first cluster => same behavior as alignBestCluster
            # firstcluster = clusters[0]

            # Get best cluster by length-normalized best score.
            #   Length normalization divides the score by the expected probability
            #   values if all peakgroups were chosen randomly (assuming equal
            #   probability between 0 and aligned_fdr_cutoff, the expected value
            #   for a random peakgroup is "aligned_fdr_cutoff/2") and thus the
            #   expected random value of n peakgroups would be (aligned_fdr_cutoff/2)^n
            clusters.sort(lambda x,y: 
                          cmp(x.getTotalScore()/((self._aligned_fdr_cutoff/2)**len(x.peakgroups)),
                          y.getTotalScore()/((self._aligned_fdr_cutoff/2)**len(y.peakgroups))) )
            # bestcluster = cluster[0]
            for i,c in enumerate(clusters): 
                if self.verbose:
                    print " - Cluster", i, "with score", c.getTotalScore(), "at", \
                      c.getMedianRT(), "+/-", c.getRTstd() , "(norm_score %s)" %\
                      (float(c.getTotalScore())/((self._aligned_fdr_cutoff/2)**len(c.peakgroups)))
                for pg in c.peakgroups:
                    pg.setClusterID(i+1)
                    if self.verbose:
                        print "   = Have member", pg.print_out()

    def _findAllPGForSeed(self, tree, tr_data, m, seed, already_seen):
        """Align peakgroups against the given seed.

        Using the given seed, the algorithm will traverse the MST tree and add
        the best matching peakgroup of that node to the result.

        Args:
            tree(list(tuple)): a minimum spanning tree (MST) represented as
                list of edges (for example [('0', '1'), ('1', '2')] ). Node names
                need to correspond to run ids.
            tr_data(format.TransformationCollection.LightTransformationData):
                structure to hold binary transformations between two different
                retention time spaces
            m(Multipeptide): one multipeptides on which the alignment should be performed
            seed(PeakGroupBase): one peakgroup chosen as the seed
            already_seen(dict): list of peakgroups already aligned (e.g. in a
                previous cluster) and which should be ignored

        Returns:
            list(PeakGroupBase): List of peakgroups belonging to this cluster
        """

        seed_rt = seed.get_normalized_retentiontime()

        # Keep track of which nodes we have already visited in the graph
        # (also storing the rt at which we found the signal in this run).
        rt_map = { seed.peptide.run.get_id() : seed_rt } 
        visited = { seed.peptide.run.get_id() : seed } 

        while len(visited.keys()) < m.get_nr_runs():
            for e1, e2 in tree:
                if e1 in visited.keys() and not e2 in visited.keys():
                    # print "try to align", e2, "from", e1
                    newPG, rt = self._findBestPG(m, e1, e2, tr_data, rt_map[e1], already_seen)
                    rt_map[e2] = rt
                    visited[e2] = newPG
                if e2 in visited.keys() and not e1 in visited.keys():
                    # print "try to align", e1, "from", e2
                    newPG, rt = self._findBestPG(m, e2, e1, tr_data, rt_map[e2], already_seen)
                    rt_map[e1] = rt
                    visited[e1] = newPG
        return [pg for pg in visited.values() if pg is not None]

    def _findBestPG(self, m,  source, target, tr_data, source_rt, already_seen):
        """Find (best) matching peakgroup in "target" which matches to the source_rt RT.

        Args:
            m(Multipeptide): one multipeptides on which the alignment should be performed
            source(string): id of the source run (where RT is known)
            target(string): id of the target run (in which the best peakgroup should be selected)
            tr_data(format.TransformationCollection.LightTransformationData):
                structure to hold binary transformations between two different
                retention time spaces
            source_rt(float): retention time of the correct peakgroup in the source run
            already_seen(dict): list of peakgroups already aligned (e.g. in a
                previous cluster) and which should be ignored

        Returns:
            list(PeakGroupBase): List of peakgroups belonging to this cluster
        """
        # Get expected RT (transformation of source into target domain)
        expected_rt = tr_data.getTrafo(source, target).predict([source_rt])[0]

        # If there is no peptide present in the target run, we simply return
        # the expected retention time.
        if not m.has_peptide(target):
            return None, expected_rt

        max_rt_diff = self._max_rt_diff
        if self._stdev_max_rt_per_run is not None:
            max_rt_diff = self._stdev_max_rt_per_run * tr_data.getStdev(source, target)
            max_rt_diff = max(self._max_rt_diff, max_rt_diff)

        # Select matching peakgroups from the target run (within the user-defined maximal rt deviation)
        target_p = m.get_peptide(target)
        matching_peakgroups = [pg_ for pg_ in target_p.get_all_peakgroups() 
            if (abs(float(pg_.get_normalized_retentiontime()) - float(expected_rt)) < max_rt_diff) and
                pg_.get_fdr_score() < self._aligned_fdr_cutoff and 
                pg_.get_feature_id() + pg_.peptide.get_id() not in already_seen]

        # If there are no peak groups present in the target run, we simply
        # return the expected retention time.
        if len(matching_peakgroups) == 0:
            return None, expected_rt

        # Select best scoring peakgroup among those in the matching RT window
        bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))

        # Printing for debug mode
        if self.verbose:
            closestPG = min(matching_peakgroups, key=lambda x: abs(float(x.get_normalized_retentiontime()) - expected_rt))
            print "  closest:", closestPG.print_out(), "diff", abs(closestPG.get_normalized_retentiontime() - expected_rt)
            print "  bestScoring:", bestScoringPG.print_out(), "diff", abs(bestScoringPG.get_normalized_retentiontime() - expected_rt)
            print

        # Decide which retention time to return:
        #  - the threading one based on the alignment
        #  - the one of the best peakgroup
        if self._correctRT_using_pg:
            return bestScoringPG, bestScoringPG.get_normalized_retentiontime()
        else:
            return bestScoringPG, expected_rt

