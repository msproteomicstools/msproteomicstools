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
def static_cy_findAllPGForSeed(tree, tr_data, multip, CyPeakgroupWrapperOnly seed, 
        dict _already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
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

    ## print ("best seed:", seed, hex(id(seed)))
    ## print ("memaddr: {0:x}".format(<long> seed.inst) )

    cdef libcpp_map[libcpp_string, int] already_seen

    # Keep track of which nodes we have already visited in the graph
    # (also storing the rt at which we found the signal in this run).
    ### cdef dict rt_map = { seed.getPeptide().getRunId() : seed_rt }
    ### cdef dict visited = { seed.getPeptide().getRunId() : seed } 

    cdef libcpp_map[libcpp_string, c_peakgroup*] c_visited
    cdef libcpp_map[libcpp_string, double] c_rt_map
    # cdef libcpp_vector[c_peakgroup].iterator pg_it

    c_visited[ seed.getPeptide().getRunId() ] = seed.inst
    c_rt_map[ seed.getPeptide().getRunId() ] = seed_rt

    cdef libcpp_pair[double, c_peakgroup*] retval

    cdef c_peakgroup * newPG
    cdef double rt
    cdef int nr_runs = multip.get_nr_runs()

    # Tree is a list of tuple( str, str)
    cdef libcpp_string e1
    cdef libcpp_string e2
    cdef libcpp_vector[ libcpp_pair[ libcpp_string, libcpp_string] ] c_tree = tree
    cdef libcpp_vector[ libcpp_pair[ libcpp_string, libcpp_string] ].iterator tree_it
    ## for e1, e2 in tree:
    ##     s1 = <bytes>e1
    ##     s2 = <bytes>e2
    ##     c_pair = 


    ## while len(visited.keys()) < nr_runs:
    while c_visited.size() < nr_runs:
        # print ("working on tree", c_visited.size(), " with totl size", nr_runs)
        tree_it = c_tree.begin()
        # for e1, e2 in tree:
        while tree_it != c_tree.end():
            e1 = deref(tree_it).first
            e2 = deref(tree_it).second
      
            #if e1 in visited.keys() and not e2 in visited.keys():
            if c_visited.find(e1) != c_visited.end() and c_visited.find(e2) == c_visited.end():
                if verbose:
                    print("  try to align", e2, "from already known node", e1)
                newPG = NULL
                retval = static_cy_findBestPG(multip, e1, e2, tr_data, c_rt_map[e1], 
                                    already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
                                    max_rt_diff, stdev_max_rt_per_run, use_local_stdev, rt,
                                    verbose)
                        
                newPG = retval.second
                rt = retval.first
                c_visited[ e2 ] = newPG
                c_rt_map[ e2 ] = rt
                # rt_map[e2] = rt
                # visited[e2] = newPG
            # if e2 in visited.keys() and not e1 in visited.keys():
            if c_visited.find(e2) != c_visited.end() and c_visited.find(e1) == c_visited.end():
                if verbose:
                    print( "  try to align", e1, "from", e2)
                newPG = NULL
                retval = static_cy_findBestPG(multip, e2, e1, tr_data, c_rt_map[e2], 
                                    already_seen, aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
                                    max_rt_diff, stdev_max_rt_per_run, use_local_stdev, rt,
                                    verbose)
                newPG = retval.second
                rt = retval.first
                # rt_map[e1] = rt
                # visited[e1] = newPG
                c_rt_map[ e1 ] = rt
                c_visited[ e1 ] = newPG

            inc(tree_it)

    ## print ("---- done with tree")
    # Now in each run at most one (zero or one) peakgroup got selected for
    # the current peptide label group. This means that for each run, either
    # heavy or light was selected (but not both) and now we should align
    # the other isotopic channels as well.
    if verbose:
        print( "Re-align isotopic channels:")

    isotopically_added_pg = []
    cdef libcpp_vector[c_peakgroup*] isotopically_added_pg_c
    cdef libcpp_map[libcpp_string, c_peakgroup*].iterator pg_it = c_visited.begin()
    cdef c_peakgroup* pg_
    cdef c_precursor* ref_peptide
    cdef CyPrecursorWrapperOnly pep
    # for pg in visited.values():
    while pg_it != c_visited.end():
        pg_ = deref(pg_it).second
        if pg_ != NULL:

            # Iterate through all sibling peptides (same peptide but
            # different isotopic composition) and pick a peak using the
            # already aligned channel as a reference.
            ref_peptide = deref(pg_).getPeptide()
            run_id = deref(ref_peptide).getRunId()

            # Iterate over CyPrecursorGroup
            # print("Try and iterate here:")
            for pep_ in multip.getPrecursorGroup(run_id):
                # print("pep1 got it, will try to use", type(pep_))
                pep = pep_
                # print (" iter here on pep", pep.get_id_c(), pep_, hex(id(pep_)), pep, hex(id(pep)),  
                #         "memaddr: {0:x}".format(<long> pep.inst))
                #       "memaddr: {0:x}".format(<long> pep_.inst) )
                # print("pep1.1 got it")
                # print("pep1.1 got it", <bytes>pep.get_id_c())
                if deref(ref_peptide).get_id() != pep.get_id_c():
                    # print("pep1.2 enter ")
                    if verbose: 
                        pass
                        # print("  Using reference %s at RT %s to align peptide %s." % (ref_peptide, pg.get_normalized_retentiontime(), pep))
                    newPG = NULL
                    # print("pep1.3 enter ")
                    # TODO: pep is of type CyPrecursorWrapperOnly -- but we expect CyPrecursorGroup 
                    retval = static_cy_findBestPGFromTemplate(deref(pg_).normalized_retentiontime, pep, max_rt_diff_isotope, already_seen,
                                                          aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, rt, verbose)
                                                          # aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, address(newPG), verbose)

                    newPG = retval.second
                    rt = retval.first
                    ### # isotopically_added_pg.append(newPG)
                    ### # print("pep1.4 enter ")
                    ### if newPG == NULL:
                    ###     print("its null here")
                    ###     print("put in pep ", type(pep))
                    ### else:
                    ###     print ("not null here")
                    ###     print ("marking pg", deref(newPG).internal_id_)

                    ## that should be nr 25
                    isotopically_added_pg_c.push_back(newPG)

                else:
                    pass
                    # print("pep1.3 dontenter ")

                # print("pep_end1 got it")
        inc(pg_it)

    if verbose:
        print("Done with re-alignment of isotopic channels")

    res = []
    cdef libcpp_vector[c_peakgroup*].iterator pg_it2 = isotopically_added_pg_c.begin()
    # TODO : something is notright here: segfault!!! 
    if False:
        pg_it = c_visited.begin()
        while pg_it != c_visited.end():
            pg_ = deref(pg_it).second
            if pg_ != NULL:

                # def __init__(self, bytes this_id, run):
                result_pep = CyPrecursorWrapperOnly(None, None)
                result_pep.inst = deref(pg_).getPeptide()

                if result_pep.inst == NULL:
                    print ("is null!")
                else:
                    print ("is not null!")

                print ("size",  deref(result_pep.inst).peakgroups.size())

                result = CyPeakgroupWrapperOnly()
                result.inst = pg_
                result.peptide = result_pep
                res.append(result)

            inc(pg_it)

        pg_it2 = isotopically_added_pg_c.begin()
        while pg_it2 != isotopically_added_pg_c.end():
            pg_ = deref(pg_it2)
            if pg_ != NULL:
                result_pep = CyPrecursorWrapperOnly(None, None)
                result_pep.inst = deref(pg_).getPeptide()

                result = CyPeakgroupWrapperOnly()
                result.inst = pg_
                result.peptide = result_pep
                res.append(result)

            inc(pg_it2)
    else:
        pg_it = c_visited.begin()
        while pg_it != c_visited.end():
            pg_ = deref(pg_it).second
            if pg_ != NULL:
                deref(pg_).cluster_id_ = 1
                ## print (" a) mark peak at : ", "memaddr: {0:x}".format(<long> pg_))
            inc(pg_it)

        pg_it2 = isotopically_added_pg_c.begin()
        while pg_it2 != isotopically_added_pg_c.end():
            pg_ = deref(pg_it2)
            if pg_ != NULL:
                deref(pg_).cluster_id_ = 1
                ref_peptide = deref(pg_).getPeptide()
                # print (" b) mark peak at : ", "memaddr: {0:x}".format(<long> pg_))

                ## print (" corresponding pep : ", "memaddr: {0:x}".format(<long> ref_peptide))
                ## ### TODO: now check that the ones from this pep are really makred!! ... 
                ## print (" corresponding pep : ", ref_peptide.curr_id_)
                ### ref_peptide.printAddresses()
                ## ref_peptide.sequence_ = "AATEST"

                ## result_pep = CyPrecursorWrapperOnly(None, None, False)
                ## result_pep.inst = deref(pg_).getPeptide()
                ## result_pep.printAddresses()

            inc(pg_it2)

    ##print("===========================================================================")
    ##print("===========================================================================")

    # return [pg for pg in list(visited.values()) + isotopically_added_pg if pg is not None]
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef libcpp_pair[double, c_peakgroup*] static_cy_findBestPG(multip, libcpp_string source, libcpp_string target, 
        CyLightTransformationData tr_data, double source_rt, 
        libcpp_map[libcpp_string, int] already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        double max_rt_diff, double stdev_max_rt_per_run, bool use_local_stdev, double &rt,
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
    ## print ("static_cy_findBestPG started")
    # Get expected RT (transformation of source into target domain)
    cdef CyLinearInterpolateWrapper cytrafo = tr_data.getTrafoCy(source, target)
    cdef libcpp_pair[double, c_peakgroup*] retval
    retval.second = NULL
    cdef double expected_rt = cytrafo.predict_cy(source_rt)
    # expected_rt = tr_data.getTrafo(source, target).predict([source_rt])[0]

    ## if True:
    ##     print(expected_rt, expected_rt_x)
    ## if verbose:
    ##     print("  Expected RT", expected_rt, " (source RT)", source_rt )
    ##     print("  --- and back again :::  ", tr_data.getTrafo(target, source).predict([expected_rt])[0] )

    # If there is no peptide present in the target run, we simply return
    # the expected retention time.
    cdef bool has_gr = multip.hasPrecursorGroup(target)
    if not has_gr:
        # return None, expected_rt
        retval.first = expected_rt
        return retval

    if stdev_max_rt_per_run is not None:
        max_rt_diff = stdev_max_rt_per_run * tr_data.getStdevCy(source, target)

        # Whether to use the standard deviation in this local region of the
        # chromatogram (if available)
        if use_local_stdev:
            max_rt_diff = stdev_max_rt_per_run * tr_data.getTrafo(source, target).last_dispersion

        max_rt_diff = max(max_rt_diff, max_rt_diff)

    ## if verbose:
    ##     print("  Used rt diff:", max_rt_diff)

    cdef CyPrecursorGroup target_peptide = multip.getPrecursorGroup(target)
    return static_cy_findBestPGFromTemplate(expected_rt, target_peptide, max_rt_diff, already_seen, 
             aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg, rt, verbose)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef libcpp_pair[double, c_peakgroup*] static_cy_findBestPGFromTemplate(double expected_rt, target_peptide, double max_rt_diff,
        libcpp_map[libcpp_string, int] already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg, double &rt,
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

    """
    # print ("static_cy_findBestPGFromTemplate started")
    # Select matching peakgroups from the target run (within the user-defined maximal rt deviation)
    cdef libcpp_pair[double, c_peakgroup*] retval
    retval.second = NULL

    cdef libcpp_vector[c_peakgroup].iterator pg_it
    cdef libcpp_vector[c_peakgroup*] matching_pg 
    cdef CyPrecursorWrapperOnly cypr
    cdef c_precursor * pep
    cdef libcpp_string id_string
    ## print ("with target peptide", target_peptide, type(target_peptide))
    if True:
        # print ("get going!", type(target_peptide))
        for precursor in target_peptide.getAllPrecursors():
            # print ("iter! going!")
            cypr = precursor
            # print ("copied ??  !")
            pg_it = cypr.get_all_peakgroups_cy_begin()
            # print ("before: getPeptidePrecursor !")
            pep = cypr.getPeptidePrecursor()
            # print ("start second iter!")
            while pg_it != cypr.get_all_peakgroups_cy_end():
                # if (abs(float(pg_.get_normalized_retentiontime()) - float(expected_rt)) < max_rt_diff) and
                if abs(deref(pg_it).normalized_retentiontime - expected_rt) < max_rt_diff:
                    # pg_.get_fdr_score() < aligned_fdr_cutoff and 
                    if deref(pg_it).fdr_score < aligned_fdr_cutoff:
                        # pg_.get_feature_id() + pg_.getPeptide().get_id() not in already_seen]
                        id_string = deref(pg_it).internal_id_ + deref(pep).curr_id_
                        ## if id_string not in already_seen:
                        if already_seen.find(id_string) == already_seen.end():
                            matching_pg.push_back( address(deref(pg_it) ))
                        # pg_.get_feature_id() + pg_.getPeptide().get_id() not in already_seen]
                inc(pg_it)

    # print (" matching pg length", matching_pg.size() )

    # If there are no peak groups present in the target run, we simply
    # return the expected retention time.
    ### if len(matching_peakgroups) == 0:
    if matching_pg.empty():
        # return None, expected_rt
        retval.first = expected_rt
        return retval

    # Select best scoring peakgroup among those in the matching RT window
    ### bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))
    cdef libcpp_vector[ c_peakgroup* ].iterator iter_m = matching_pg.begin()
    cdef libcpp_vector[ c_peakgroup* ].iterator bestScoring_pg = matching_pg.begin()
    cdef double best_fdr = deref(iter_m).fdr_score
    while iter_m != matching_pg.end():
        if (deref(iter_m).fdr_score <= best_fdr):
            best_fdr = deref(iter_m).fdr_score
            bestScoring_pg = iter_m
        inc(iter_m)

    # print (" best scoring ", deref(bestScoring_pg).internal_id_)

    #### # Printing for debug mode

    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._aligned_fdr_cutoff]) > 1:
    ###     self.nr_multiple_align += 1
    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._fdr_cutoff]) > 1:
    ###     self.nr_ambiguous += 1

    # Decide which retention time to return:
    #  - the threading one based on the alignment
    #  - the one of the best peakgroup

    ## if correctRT_using_pg:
    ##     return deref(bestScoring_pg), deref(bestScoring_pg).normalized_retentiontime
    ## else:
    ##     return deref(bestScoring_pg), expected_rt
    # deref(best_pg) = 
    retval.second = deref(bestScoring_pg)

    if correctRT_using_pg:
        # rt = deref(bestScoring_pg).normalized_retentiontime
        retval.first = deref(bestScoring_pg).normalized_retentiontime
        return retval
    else:
        retval.first = expected_rt
        return retval

