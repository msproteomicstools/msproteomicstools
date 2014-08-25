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
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree
import msproteomicstoolslib.math.Smoothing as smoothing

def getMinimumSpanningTree(exp, multipeptides, initial_alignment_cutoff):

    spl_aligner = SplineAligner(initial_alignment_cutoff)
    dist_matrix = numpy.zeros(shape=(len(exp.runs),len(exp.runs)))
    for i in range(len(exp.runs)):
        for j in range(len(exp.runs)):
            if i == j:
                dist_matrix[i,j] = 0
                continue 

            idata, jdata = spl_aligner._getRTData(exp.runs[i], exp.runs[j], multipeptides)

            # Get linear alignment
            smlin = smoothing.SmoothingLinear()
            smlin.initialize(idata, jdata)
            idata_lin_aligned = smlin.predict(idata)

            # Use stdev to estimate distance between two runs
            stdev_lin = numpy.std(numpy.array(jdata) - numpy.array(idata_lin_aligned))
            dist_matrix[i,j] = stdev_lin

    return MinimumSpanningTree(dist_matrix)

class TreeConsensusAlignment():

    def __init__(self, max_rt_diff, aligned_fdr_cutoff):
        self._max_rt_diff = max_rt_diff
        self._aligned_fdr_cutoff = aligned_fdr_cutoff

    def alignBestCluster(self, exp, multipeptides, tree, tr_data):

        verb = False
        for m in multipeptides:

            # Find the overall best peptide
            best = m.find_best_peptide_pg()
            if verb: 
                print "00000000000000000000000000000000000 new peptide (cluster)", m.get_peptides()[0].get_id()
                print " Best", best.print_out(), "from run", best.peptide.run.get_id()

            # Use this peptide to generate a cluster
            for pg_ in self._findPGCluster(tree, tr_data, m, best, {}):
                pg_.select_this_peakgroup()

    def alignAllCluster(self, exp, multipeptides, tree, tr_data):

        from msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm import Cluster
        for m in multipeptides:

            verb = False
            last_cluster = []
            already_seen = set([])
            stillLeft = [b for a in m.get_peptides() for b in a.get_all_peakgroups() if b.get_feature_id() + b.peptide.get_id() not in already_seen and b.get_fdr_score() < 0.01]
            clusters = []
            while len(stillLeft) > 0:
                best = min(stillLeft, key=lambda x: float(x.get_fdr_score()))
                last_cluster = self._findPGCluster(tree, tr_data, m, best, already_seen)
                already_seen.update( set([ b.get_feature_id() + b.peptide.get_id() for b in last_cluster if b is not None]) )
                stillLeft = [ b for a in m.get_peptides() for b in a.get_all_peakgroups() if b.get_feature_id() + b.peptide.get_id() not in already_seen and b.get_fdr_score() < 0.01]
                clusters.append(Cluster(last_cluster))

            # select the first cluster => same behavior as alignBestCluster
            # firstcluster = clusters[0]

            # Get best cluster by length-normalized best score.
            #   Length normalization divides the score by the expected probability
            #   values if all peakgroups were chosen randomly (assuming equal
            #   probability between 0 and aligned_fdr_cutoff, the expected value
            #   for a random peakgroup is "aligned_fdr_cutoff/2") and thus the
            #   expected random value of n peakgroups would be (aligned_fdr_cutoff/2)^n
            cluster_order = clusters.sort(lambda x,y: 
                                          cmp(x.getTotalScore()/(((self._aligned_fdr_cutoff/2)**len(x.peakgroups))),
                                              y.getTotalScore()/(((self._aligned_fdr_cutoff/2)**len(y.peakgroups)))) )
            # bestcluster = cluster[0]
            for i,c in enumerate(clusters): 
                if False:
                    print " - Cluster", i, "with score", c.getTotalScore(), "at", \
                      c.getMedianRT(), "+/-", c.getRTstd() , "(norm_score %s)" %\
                      (float(c.getTotalScore())/((self._aligned_fdr_cutoff/2)**len(c.peakgroups)))
                for pg in c.peakgroups:
                    pg.setClusterID(i+1)
                    if False:
                        print "   = Have member", pg.print_out()

    def _findPGCluster(self, tree, tr_data, m, seed, already_seen):
        seed_rt = seed.get_normalized_retentiontime()

        # Keep track of which nodes we have already visited in the graph
        # (also storing the rt at which we found the signal in this run).
        rt_map = { seed.peptide.run.get_id() : seed_rt } 
        visited = { seed.peptide.run.get_id() : seed } 

        while len(visited.keys()) != m.get_nr_runs():
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

        # Get expected RT (transformation of source into target domain)
        expected_rt = tr_data.getTrafo(source, target).predict([source_rt])[0]

        # If there is no peptide present in the target run, we simply return
        # the expected retention time.
        if not m.has_peptide(target):
            return None, expected_rt

        # Select matching peakgroups from the target run (within the user-defined maximal rt deviation)
        target_p = m.get_peptide(target)
        matching_peakgroups = [pg_ for pg_ in target_p.get_all_peakgroups() 
            if (abs(float(pg_.get_normalized_retentiontime()) - float(expected_rt)) < self._max_rt_diff) and
                pg_.get_fdr_score() < self._aligned_fdr_cutoff and 
                pg_.get_feature_id() + pg_.peptide.get_id() not in already_seen ]

        # If there are no peak groups present in the target run, we simply
        # return the expected retention time.
        if len(matching_peakgroups) == 0:
            return None, expected_rt

        # Select best scoring peakgroup among those in the matching RT window
        bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))
        # closestPG = min(matching_peakgroups, key=lambda x: abs(float(x.get_normalized_retentiontime()) - expected_rt))
        # print "  closest:", closestPG.print_out(), "diff", abs(closestPG.get_normalized_retentiontime() - expected_rt)
        # print "  bestScoring:", bestScoringPG.print_out(), "diff", abs(bestScoringPG.get_normalized_retentiontime() - expected_rt)
        # print
        return bestScoringPG, expected_rt

