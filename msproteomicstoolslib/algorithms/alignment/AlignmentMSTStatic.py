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

from __future__ import print_function
import numpy
import math
import scipy.stats
from msproteomicstoolslib.format.TransformationCollection import LightTransformationData
from msproteomicstoolslib.data_structures.PrecursorGroup import PrecursorGroup
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
import msproteomicstoolslib.data_structures.PeakGroup
import msproteomicstoolslib.math.Smoothing as smoothing

def static_findAllPGForSeed(tree, tr_data, m, seed, 
        already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
        max_rt_diff, stdev_max_rt_per_run, use_local_stdev, max_rt_diff_isotope,
        verbose):
    """Align peakgroups against the given seed.

    Using the given seed, the algorithm will traverse the MST tree and add
    the best matching peakgroup of that node to the result.

    Args:
        tree(list(tuple)): a minimum spanning tree (MST) represented as
            list of edges (for example [('0', '1'), ('1', '2')] ). Node names
            need to correspond to run ids.
        tr_data(:class:`.LightTransformationData`): structure to hold
            binary transformations between two different retention time spaces
        m(Multipeptide): one multipeptides on which the alignment should be performed
        seed(PeakGroupBase): one peakgroup chosen as the seed
        already_seen(dict): list of peakgroups already aligned (e.g. in a
            previous cluster) and which should be ignored

    Returns:
        list(PeakGroupBase): List of peakgroups belonging to this cluster
    """
    ### print ("def static_findAllPGForSeed(tree, tr_data, m, seed, ", test_optimized(5, 6) )
    if verbose: 
        print("111111111111111111111111 findAllPGForSeed started" )
        print("  Seed", seed.print_out(), "from run", seed.getPeptide().getRunId() )

    seed_rt = seed.get_normalized_retentiontime()

    # Keep track of which nodes we have already visited in the graph
    # (also storing the rt at which we found the signal in this run).
    rt_map = { seed.getPeptide().getRunId() : seed_rt } 
    visited = { seed.getPeptide().getRunId() : seed } 

    while len(visited.keys()) < m.get_nr_runs():
        for e1, e2 in tree:
            if e1 in visited.keys() and not e2 in visited.keys():
                if verbose:
                    print("  try to align", e2, "from already known node", e1)
                newPG, rt = static_findBestPG(m, e1, e2, tr_data, rt_map[e1], 
                                    already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
                                    max_rt_diff, stdev_max_rt_per_run, use_local_stdev,
                                    verbose)
                        
                rt_map[e2] = rt
                visited[e2] = newPG
            if e2 in visited.keys() and not e1 in visited.keys():
                if verbose: 
                    print( "  try to align", e1, "from", e2)
                newPG, rt = static_findBestPG(m, e2, e1, tr_data, rt_map[e2], 
                                    already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
                                    max_rt_diff, stdev_max_rt_per_run, use_local_stdev,
                                    verbose)
                rt_map[e1] = rt
                visited[e1] = newPG

    # Now in each run at most one (zero or one) peakgroup got selected for
    # the current peptide label group. This means that for each run, either
    # heavy or light was selected (but not both) and now we should align
    # the other isotopic channels as well.
    if verbose: 
        print( "Re-align isotopic channels:")

    isotopically_added_pg = []
    for pg in visited.values():
        if pg is not None:

            # Iterate through all sibling peptides (same peptide but
            # different isotopic composition) and pick a peak using the
            # already aligned channel as a reference.
            ref_peptide = pg.getPeptide()
            run_id = ref_peptide.getRunId()
            for pep in m.getPrecursorGroup(run_id):
                if ref_peptide != pep:
                    if verbose: 
                        print("  Using reference %s at RT %s to align peptide %s." % (ref_peptide, pg.get_normalized_retentiontime(), pep))
                    newPG, rt = static_findBestPGFromTemplate(pg.get_normalized_retentiontime(), pep, max_rt_diff_isotope, already_seen,
                                                              aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, verbose)
                    isotopically_added_pg.append(newPG)

    if verbose: 
        print("Done with re-alignment of isotopic channels")

    return [pg for pg in list(visited.values()) + isotopically_added_pg if pg is not None]

def static_findBestPG(m, source, target, tr_data, source_rt, 
        already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
        max_rt_diff, stdev_max_rt_per_run, use_local_stdev,
        verbose):
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
    if verbose:
        print("  Expected RT", expected_rt, " (source RT)", source_rt )
        print("  --- and back again :::  ", tr_data.getTrafo(target, source).predict([expected_rt])[0] )

    # If there is no peptide present in the target run, we simply return
    # the expected retention time.
    if not m.hasPrecursorGroup(target):
        return None, expected_rt

    if stdev_max_rt_per_run is not None:
        max_rt_diff = stdev_max_rt_per_run * tr_data.getStdev(source, target)

        # Whether to use the standard deviation in this local region of the
        # chromatogram (if available)
        if use_local_stdev:
            max_rt_diff = stdev_max_rt_per_run * tr_data.getTrafo(source, target).last_dispersion

        max_rt_diff = max(max_rt_diff, max_rt_diff)

    if verbose:
        print("  Used rt diff:", max_rt_diff)

    ### print (type( source ) )
    ### print (type( m ) )
    ### print (type( m.getPrecursorGroup(target) ) )
    return static_findBestPGFromTemplate(expected_rt, m.getPrecursorGroup(target), max_rt_diff, already_seen, 
             aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, verbose)

def static_findBestPGFromTemplate(expected_rt, target_peptide, max_rt_diff,
        already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
        verbose):
    """Find (best) matching peakgroup in "target" which matches to the source_rt RT.

        Parameters
        ----------
        expected_rt : float
            Expected retention time
        target_peptide: :class:`.PrecursorGroup`
            Precursor group from the target run (contains multiple peak groups)
        max_rt_diff : float
            Maximal retention time difference (parameter)
        already_seen : dict
            list of peakgroups already aligned (e.g. in a previous cluster) and which should be ignored
        aligned_fdr_cutoff : float
        fdr_cutoff : float
        correctRT_using_pg: boolean
        verbose: boolean
    """
    # Select matching peakgroups from the target run (within the user-defined maximal rt deviation)
    matching_peakgroups = [pg_ for pg_ in target_peptide.getAllPeakgroups() 
        if (abs(float(pg_.get_normalized_retentiontime()) - float(expected_rt)) < max_rt_diff) and
            pg_.get_fdr_score() < aligned_fdr_cutoff and 
            pg_.get_feature_id() + pg_.getPeptide().get_id() not in already_seen]

    # If there are no peak groups present in the target run, we simply
    # return the expected retention time.
    if len(matching_peakgroups) == 0:
        return None, expected_rt

    # Select best scoring peakgroup among those in the matching RT window
    bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))

    # Printing for debug mode
    if verbose:
        closestPG = min(matching_peakgroups, key=lambda x: abs(float(x.get_normalized_retentiontime()) - expected_rt))
        print("    closest:", closestPG.print_out(), "diff", abs(closestPG.get_normalized_retentiontime() - expected_rt) )
        print("    bestScoring:", bestScoringPG.print_out(), "diff", abs(bestScoringPG.get_normalized_retentiontime() - expected_rt) )
        print()

    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._aligned_fdr_cutoff]) > 1:
    ###     self.nr_multiple_align += 1
    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._fdr_cutoff]) > 1:
    ###     self.nr_ambiguous += 1

    # Decide which retention time to return:
    #  - the threading one based on the alignment
    #  - the one of the best peakgroup
    if correctRT_using_pg:
        return bestScoringPG, bestScoringPG.get_normalized_retentiontime()
    else:
        return bestScoringPG, expected_rt


