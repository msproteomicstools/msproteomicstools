# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np

from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from libcpp cimport bool

import numpy as np
import operator

@cython.boundscheck(False)
@cython.wraparound(False)
def static_findAllPGForSeed(tree, tr_data, multip, CyPeakgroupWrapperOnly seed, 
        dict already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        double max_rt_diff, double stdev_max_rt_per_run, bool use_local_stdev, double max_rt_diff_isotope,
        bool verbose):
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

    seed_rt = seed.get_normalized_retentiontime_cy()

    # Keep track of which nodes we have already visited in the graph
    # (also storing the rt at which we found the signal in this run).
    cdef dict rt_map = { seed.getPeptide().getRunId() : seed_rt }
    cdef dict visited = { seed.getPeptide().getRunId() : seed } 

    cdef int nr_runs = multip.get_nr_runs()
    while len(visited.keys()) < nr_runs:
        for e1, e2 in tree:
            if e1 in visited.keys() and not e2 in visited.keys():
                if verbose:
                    print("  try to align", e2, "from already known node", e1)
                newPG, rt = static_findBestPG(multip, e1, e2, tr_data, rt_map[e1], 
                                    already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
                                    max_rt_diff, stdev_max_rt_per_run, use_local_stdev,
                                    verbose)
                        
                rt_map[e2] = rt
                visited[e2] = newPG
            if e2 in visited.keys() and not e1 in visited.keys():
                if verbose: 
                    print( "  try to align", e1, "from", e2)
                newPG, rt = static_findBestPG(multip, e2, e1, tr_data, rt_map[e2], 
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

            for pep in multip.getPrecursorGroup(run_id):
                if ref_peptide.get_id() != pep.get_id():
                    if verbose: 
                        print("  Using reference %s at RT %s to align peptide %s." % (ref_peptide, pg.get_normalized_retentiontime(), pep))
                    newPG, rt = static_findBestPGFromTemplate(pg.get_normalized_retentiontime(), pep, max_rt_diff_isotope, already_seen,
                                                                aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, verbose)
                    isotopically_added_pg.append(newPG)
                    # print ("marking pg", newPG, newPG.get_feature_id())

    if verbose: 
        print("Done with re-alignment of isotopic channels")

    return [pg for pg in list(visited.values()) + isotopically_added_pg if pg is not None]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef static_findBestPG(multip, source, target, CyLightTransformationData tr_data, double source_rt, 
        dict already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        double max_rt_diff, double stdev_max_rt_per_run, bool use_local_stdev,
        bool verbose):
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
    cdef CyLinearInterpolateWrapper cytrafo = tr_data.getTrafoCy(source, target)
    cdef double expected_rt = cytrafo.predict_cy(source_rt)
    # expected_rt = tr_data.getTrafo(source, target).predict([source_rt])[0]

    ## if True:
    ##     print(expected_rt, expected_rt_x)

    if verbose:
        print("  Expected RT", expected_rt, " (source RT)", source_rt )
        print("  --- and back again :::  ", tr_data.getTrafo(target, source).predict([expected_rt])[0] )

    # If there is no peptide present in the target run, we simply return
    # the expected retention time.
    cdef bool has_gr = multip.hasPrecursorGroup(target)
    if not has_gr:
        return None, expected_rt

    if stdev_max_rt_per_run is not None:
        max_rt_diff = stdev_max_rt_per_run * tr_data.getStdevCy(source, target)

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
    return static_findBestPGFromTemplate(expected_rt, multip.getPrecursorGroup(target), max_rt_diff, already_seen, 
             aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, verbose)
    ## return static_cy_findBestPGFromTemplate(expected_rt, m.getPrecursorGroup(target), max_rt_diff, already_seen, 
    ##         aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, verbose)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef static_findBestPGFromTemplate(double expected_rt, target_peptide, double max_rt_diff,
        dict already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        bool verbose):
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

    # print (" matching pg length", len(matching_peakgroups) )

    # Select best scoring peakgroup among those in the matching RT window
    bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))
    # print (" best scoring ", bestScoringPG)

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



