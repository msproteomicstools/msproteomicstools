# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np

from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.unordered_set cimport unordered_set as libcpp_unordered_set
from libcpp.pair cimport pair as libcpp_pair
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from libcpp cimport bool

## for stand-alone compile
## include "LightTransformationData.pyx"
## include "PeakgroupWrapper.pyx"
## include "PrecursorWrapper.pyx"
## include "PrecursorGroup.pyx"

cdef cppclass mst_settings:

    double aligned_fdr_cutoff 
    double fdr_cutoff
    bool correctRT_using_pg
    double max_rt_diff
    bool use_local_stdev
    double max_rt_diff_isotope
    bool verbose
    double stdev_max_rt_per_run
    bool use_stdev_max_rt_per_run
    int nr_multiple_align 
    int nr_ambiguous 

    ###     self.nr_multiple_align += 1
    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._fdr_cutoff]) > 1:
    ###     self.nr_ambiguous += 1

cdef cppclass c_node:

    libcpp_string internal_id
    libcpp_vector[c_node*] neighbors

cdef cppclass cpp_tree:

    libcpp_map[libcpp_string, c_node*] nodes

cdef visit_tree_simple(cpp_tree tree, c_node * current, libcpp_unordered_set[libcpp_string] visited):

    # mark as visited
    visited.insert( deref(current).internal_id )
    cdef c_node* node1
    cdef c_node* node2
    cdef libcpp_string e1
    cdef libcpp_vector[c_node*].iterator n_it = deref(current).neighbors.begin()
    while (n_it != deref(current).neighbors.end()):
        if visited.find( deref(deref(n_it)).internal_id) != visited.end():
            # have visited this neighbor already
            inc(n_it)
        else:
            # have not visited this edge (current, neighbor) yet
            node1 = current
            node2 = deref(n_it)
            # do work
            print ("visit!", deref(node1).internal_id, deref(node2).internal_id)
            # recursive call
            visit_tree_simple(tree, node2, visited)
            inc(n_it)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef cy_findAllPGForSeedTreeIter(cpp_tree * tree, c_node * current, libcpp_unordered_set[libcpp_string] * visited, tr_data, multip,
        libcpp_map[libcpp_string, int] already_seen, libcpp_map[libcpp_string, double] * c_rt_map, double current_rt,
        libcpp_map[libcpp_string, c_peakgroup*] * c_visited,
        bool verbose, mst_settings * settings):

    # mark as visited
    deref(visited).insert( deref(current).internal_id )
    
    cdef c_node* node1
    cdef c_node* node2
    cdef libcpp_string e1
    cdef libcpp_string e2
    cdef libcpp_pair[double, c_peakgroup*] retval
    cdef libcpp_vector[c_node*].iterator n_it = deref(current).neighbors.begin()
    cdef c_peakgroup * newPG
    while (n_it != deref(current).neighbors.end()):
        if deref(visited).find( deref(deref(n_it)).internal_id) != deref(visited).end():
            # have visited this neighbor already
            inc(n_it)
        else:
            # have not visited this edge (current, neighbor) yet
            node1 = current
            node2 = deref(n_it)
            e1 = node1.internal_id
            e2 = node2.internal_id

            # do work
            retval = static_cy_findBestPG(multip, e1, e2, tr_data, deref(c_rt_map)[e1], already_seen, settings)

            newPG = retval.second
            rt = retval.first # target_rt
            deref(c_visited)[ e2 ] = newPG
            deref(c_rt_map)[ e2 ] = rt
            # recursive call
            cy_findAllPGForSeedTreeIter(tree, node2, visited, tr_data, multip, already_seen, c_rt_map, current_rt, c_visited, verbose, settings)

            inc(n_it)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef cpp_tree createTree(libcpp_vector[ libcpp_pair[ libcpp_string, libcpp_string] ] c_tree):
    
    cdef cpp_tree tree

    # Tree is a list of tuple( str, str)
    cdef libcpp_vector[ libcpp_pair[ libcpp_string, libcpp_string] ].iterator tree_it
    tree_it = c_tree.begin()
    cdef c_node* node1
    cdef c_node* node2
    cdef libcpp_string e1
    cdef libcpp_string e2
    while tree_it != c_tree.end():

        e1 = deref(tree_it).first
        e2 = deref(tree_it).second
        if tree.nodes.find(e1) == tree.nodes.end():
            # create new e1 node
            node1 = new c_node()
            deref(node1).internal_id = e1
            tree.nodes[e1] = node1
        if tree.nodes.find(e2) == tree.nodes.end():
            # create new e2 node
            node2 = new c_node()
            deref(node2).internal_id = e2
            tree.nodes[e2] = node2

        node1 = tree.nodes[e1]
        node2 = tree.nodes[e2]

        # Add them both as neighbors to each other
        deref(node1).neighbors.push_back(node2)
        deref(node2).neighbors.push_back(node1)

        inc(tree_it)

    return tree

@cython.boundscheck(False)
@cython.wraparound(False)
def static_cy_alignBestCluster(multipeptides, py_tree, tr_data,
        double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        double max_rt_diff, stdev_max_rt_per_run, bool use_local_stdev, double max_rt_diff_isotope,
        bool verbose):
    """Use the MST to report the first cluster containing the best peptide (overall).

    The algorithm will go through all multipeptides and mark those
    peakgroups which it deems to belong to the best peakgroup cluster (only
    the first cluster will be reported).

    Args:
        multipeptides(list of :class:`.Multipeptide`): a list of
            multipeptides on which the alignment should be performed. After
            alignment, each peakgroup that should be quantified can be
            retrieved by calling get_selected_peakgroups() on the multipeptide.
        tree(list of tuple): a minimum spanning tree (MST) represented as
            list of edges (for example [('0', '1'), ('1', '2')] ). Node names
            need to correspond to run ids.
        tr_data(format.TransformationCollection.LightTransformationData):
            structure to hold binary transformations between two different
            retention time spaces

    Returns:
        None
    """

    cdef libcpp_vector[ libcpp_pair[ libcpp_string, libcpp_string] ] c_tree = py_tree
    cdef cpp_tree tree = createTree(c_tree)

    cdef mst_settings settings
    settings.aligned_fdr_cutoff  = aligned_fdr_cutoff 
    settings.fdr_cutoff          = fdr_cutoff         
    settings.correctRT_using_pg  = correctRT_using_pg 
    settings.max_rt_diff         = max_rt_diff        
    settings.use_local_stdev     = use_local_stdev    
    settings.use_stdev_max_rt_per_run = (stdev_max_rt_per_run is not None)
    settings.max_rt_diff_isotope = max_rt_diff_isotope
    settings.verbose             = verbose            
    settings.nr_multiple_align = 0
    settings.nr_ambiguous = 0
    if settings.use_stdev_max_rt_per_run:
        settings.stdev_max_rt_per_run = stdev_max_rt_per_run    

    for m in multipeptides:

        # Find the overall best peptide
        best = m.find_best_peptide_pg()

        if best.get_fdr_score() >= fdr_cutoff:
            continue

        pg_list = static_cy_fast_findAllPGForSeed(address(tree), tr_data, m, best, {}, 
                            aligned_fdr_cutoff, fdr_cutoff, correctRT_using_pg,
                            float(max_rt_diff), stdev_max_rt_per_run, use_local_stdev, max_rt_diff_isotope,
                            verbose, address(settings))

        for pg_ in pg_list:
            pg_.select_this_peakgroup()

    return (settings.nr_multiple_align, settings.nr_ambiguous)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef static_cy_fast_findAllPGForSeed(cpp_tree * c_tree, tr_data, multip, CyPeakgroupWrapperOnly seed, 
        dict _already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        double max_rt_diff, stdev_max_rt_per_run, bool use_local_stdev, double max_rt_diff_isotope,
        bool verbose, mst_settings * settings):
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

    cdef libcpp_map[libcpp_string, int] already_seen

    # Keep track of which nodes we have already visited in the graph
    # (also storing the rt at which we found the signal in this run).
    cdef libcpp_map[libcpp_string, c_peakgroup*] c_visited
    cdef libcpp_map[libcpp_string, double] c_rt_map

    cdef libcpp_string seed_run_id = seed.getPeptide().getRunId()
    c_visited[ seed_run_id ] = seed.inst
    c_rt_map[ seed_run_id ] = seed_rt

    cdef c_node * current = c_tree.nodes[ seed_run_id ]
    cdef libcpp_unordered_set[libcpp_string] visited

    # Main routine, recursively visiting the tree
    cy_findAllPGForSeedTreeIter(c_tree, current, address(visited), tr_data, multip,
        already_seen, address(c_rt_map), seed_rt, address(c_visited), verbose, settings)

    # Now in each run at most one (zero or one) peakgroup got selected for
    # the current peptide label group. This means that for each run, either
    # heavy or light was selected (but not both) and now we should align
    # the other isotopic channels as well.
    if verbose:
        print( "Re-align isotopic channels:")

    cdef libcpp_vector[c_peakgroup*] isotopically_added_pg_c
    cdef libcpp_map[libcpp_string, c_peakgroup*].iterator pg_it = c_visited.begin()
    cdef c_peakgroup* pg_
    cdef c_precursor* ref_peptide
    cdef CyPrecursorWrapperOnly pep

    cdef libcpp_pair[double, c_peakgroup*] retval
    cdef c_peakgroup * newPG

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
            for pep_ in multip.getPrecursorGroup(run_id):
                pep = pep_
                if deref(ref_peptide).get_id() != pep.get_id_c():
                    if verbose: 
                        print("  Using reference %s at RT %s to align peptide %s." % (ref_peptide.get_id(), deref(pg_).normalized_retentiontime, pep))
                    newPG = NULL
                    retval = static_cy_findBestPGFromTemplate(deref(pg_).normalized_retentiontime, pep, max_rt_diff_isotope, already_seen, settings)
                    newPG = retval.second
                    # rt = retval.first
                    isotopically_added_pg_c.push_back(newPG)

                else:
                    pass
        inc(pg_it)

    if verbose:
        print("Done with re-alignment of isotopic channels")

    # Now we need to assemble the peak groups back together
    # We need to construct the Python objects and make sure that memory
    # management is done by the original object and not the Python GC engine.
    # Therefore we take the pointers for constructing new Python objects but
    # take care not to take ownership.
    # return [pg for pg in list(visited.values()) + isotopically_added_pg if pg is not None]
    res = []
    cdef libcpp_vector[c_peakgroup*].iterator pg_it2 = isotopically_added_pg_c.begin()
    if True:
        pg_it = c_visited.begin()
        while pg_it != c_visited.end():
            pg_ = deref(pg_it).second
            if pg_ != NULL:

                result_pep = CyPrecursorWrapperOnly(None, None, False)
                result_pep.inst = deref(pg_).getPeptide()

                result = CyPeakgroupWrapperOnly()
                result.inst = pg_
                result.peptide = result_pep
                res.append(result)

            inc(pg_it)

        pg_it2 = isotopically_added_pg_c.begin()
        while pg_it2 != isotopically_added_pg_c.end():
            pg_ = deref(pg_it2)
            if pg_ != NULL:
                result_pep = CyPrecursorWrapperOnly(None, None, False)
                result_pep.inst = deref(pg_).getPeptide()

                result = CyPeakgroupWrapperOnly()
                result.inst = pg_
                result.peptide = result_pep
                res.append(result)

            inc(pg_it2)
    else:
        # This is the other option, just mark these as selected right away...
        pg_it = c_visited.begin()
        while pg_it != c_visited.end():
            pg_ = deref(pg_it).second
            if pg_ != NULL:
                deref(pg_).cluster_id_ = 1
            inc(pg_it)

        pg_it2 = isotopically_added_pg_c.begin()
        while pg_it2 != isotopically_added_pg_c.end():
            pg_ = deref(pg_it2)
            if pg_ != NULL:
                deref(pg_).cluster_id_ = 1
                ref_peptide = deref(pg_).getPeptide()
            inc(pg_it2)

    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def static_cy_findAllPGForSeed(tree, tr_data, multip, CyPeakgroupWrapperOnly seed, 
        dict _already_seen, double aligned_fdr_cutoff, double fdr_cutoff, bool correctRT_using_pg,
        double max_rt_diff, stdev_max_rt_per_run, bool use_local_stdev, double max_rt_diff_isotope,
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

    cdef mst_settings settings
    settings.aligned_fdr_cutoff  = aligned_fdr_cutoff 
    settings.fdr_cutoff          = fdr_cutoff         
    settings.correctRT_using_pg  = correctRT_using_pg 
    settings.max_rt_diff         = max_rt_diff        
    settings.use_local_stdev     = use_local_stdev    
    settings.use_stdev_max_rt_per_run = (stdev_max_rt_per_run is not None)
    settings.max_rt_diff_isotope = max_rt_diff_isotope
    settings.verbose             = verbose            
    if settings.use_stdev_max_rt_per_run:
        settings.stdev_max_rt_per_run = stdev_max_rt_per_run    

    cdef libcpp_map[libcpp_string, int] already_seen

    # Keep track of which nodes we have already visited in the graph
    # (also storing the rt at which we found the signal in this run).
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

    ## while len(visited.keys()) < nr_runs:
    while c_visited.size() < nr_runs:
        tree_it = c_tree.begin()
        while tree_it != c_tree.end():
            e1 = deref(tree_it).first
            e2 = deref(tree_it).second
      
            #if e1 in visited.keys() and not e2 in visited.keys():
            if c_visited.find(e1) != c_visited.end() and c_visited.find(e2) == c_visited.end():
                if verbose:
                    print("  try to align", e2, "from already known node", e1)
                newPG = NULL
                retval = static_cy_findBestPG(multip, e1, e2, tr_data, c_rt_map[e1], 
                                    already_seen, address(settings))
                        
                newPG = retval.second
                rt = retval.first
                c_visited[ e2 ] = newPG
                c_rt_map[ e2 ] = rt
            # if e2 in visited.keys() and not e1 in visited.keys():
            if c_visited.find(e2) != c_visited.end() and c_visited.find(e1) == c_visited.end():
                if verbose:
                    print( "  try to align", e1, "from", e2)
                newPG = NULL
                retval = static_cy_findBestPG(multip, e2, e1, tr_data, c_rt_map[e2], 
                                    already_seen, address(settings))
                newPG = retval.second
                rt = retval.first
                c_rt_map[ e1 ] = rt
                c_visited[ e1 ] = newPG

            inc(tree_it)

    # Now in each run at most one (zero or one) peakgroup got selected for
    # the current peptide label group. This means that for each run, either
    # heavy or light was selected (but not both) and now we should align
    # the other isotopic channels as well.
    if verbose:
        print( "Re-align isotopic channels:")

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
            for pep_ in multip.getPrecursorGroup(run_id):
                pep = pep_
                if deref(ref_peptide).get_id() != pep.get_id_c():
                    if verbose: 
                        print("  Using reference %s at RT %s to align peptide %s." % (ref_peptide.get_id(), deref(pg_).normalized_retentiontime, pep))
                    newPG = NULL
                    retval = static_cy_findBestPGFromTemplate(deref(pg_).normalized_retentiontime, pep, max_rt_diff_isotope, already_seen, address(settings))
                    newPG = retval.second
                    rt = retval.first
                    isotopically_added_pg_c.push_back(newPG)

                else:
                    pass
        inc(pg_it)

    if verbose:
        print("Done with re-alignment of isotopic channels")

    # Now we need to assemble the peak groups back together
    # We need to construct the Python objects and make sure that memory
    # management is done by the original object and not the Python GC engine.
    # Therefore we take the pointers for constructing new Python objects but
    # take care not to take ownership.
    # return [pg for pg in list(visited.values()) + isotopically_added_pg if pg is not None]
    res = []
    cdef libcpp_vector[c_peakgroup*].iterator pg_it2 = isotopically_added_pg_c.begin()
    if True:
        pg_it = c_visited.begin()
        while pg_it != c_visited.end():
            pg_ = deref(pg_it).second
            if pg_ != NULL:

                result_pep = CyPrecursorWrapperOnly(None, None, False)
                result_pep.inst = deref(pg_).getPeptide()

                result = CyPeakgroupWrapperOnly()
                result.inst = pg_
                result.peptide = result_pep
                res.append(result)

            inc(pg_it)

        pg_it2 = isotopically_added_pg_c.begin()
        while pg_it2 != isotopically_added_pg_c.end():
            pg_ = deref(pg_it2)
            if pg_ != NULL:
                result_pep = CyPrecursorWrapperOnly(None, None, False)
                result_pep.inst = deref(pg_).getPeptide()

                result = CyPeakgroupWrapperOnly()
                result.inst = pg_
                result.peptide = result_pep
                res.append(result)

            inc(pg_it2)
    else:
        # This is the other option, just mark these as selected right away...
        pg_it = c_visited.begin()
        while pg_it != c_visited.end():
            pg_ = deref(pg_it).second
            if pg_ != NULL:
                deref(pg_).cluster_id_ = 1
            inc(pg_it)

        pg_it2 = isotopically_added_pg_c.begin()
        while pg_it2 != isotopically_added_pg_c.end():
            pg_ = deref(pg_it2)
            if pg_ != NULL:
                deref(pg_).cluster_id_ = 1
                ref_peptide = deref(pg_).getPeptide()
            inc(pg_it2)

    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef libcpp_pair[double, c_peakgroup*] static_cy_findBestPG(multip, libcpp_string source, libcpp_string target, 
        CyLightTransformationData tr_data, double source_rt, 
        libcpp_map[libcpp_string, int] already_seen, mst_settings * settings):

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
    cdef libcpp_pair[double, c_peakgroup*] retval
    retval.second = NULL
    cdef double expected_rt = cytrafo.predict_cy(source_rt)

    if settings.verbose:
        print("  Expected RT", expected_rt, " (source RT)", source_rt )
        print("  --- and back again :::  ", tr_data.getTrafo(target, source).predict([expected_rt])[0] )

    # If there is no peptide present in the target run, we simply return
    # the expected retention time.
    cdef bool has_gr = multip.hasPrecursorGroup(target)
    if not has_gr:
        retval.first = expected_rt
        return retval
        # return None, expected_rt

    cdef double max_rt_diff = settings.max_rt_diff
    if settings.use_stdev_max_rt_per_run:
        max_rt_diff = settings.stdev_max_rt_per_run * tr_data.getStdevCy(source, target)

        # Whether to use the standard deviation in this local region of the
        # chromatogram (if available)
        if settings.use_local_stdev:
            max_rt_diff = settings.stdev_max_rt_per_run * tr_data.getTrafo(source, target).last_dispersion

        max_rt_diff = max(max_rt_diff, max_rt_diff)

    if settings.verbose:
        print("  Used rt diff:", max_rt_diff)

    cdef CyPrecursorGroup target_peptide = multip.getPrecursorGroup(target)
    return static_cy_findBestPGFromTemplate(expected_rt, target_peptide, max_rt_diff, already_seen, settings)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef libcpp_pair[double, c_peakgroup*] static_cy_findBestPGFromTemplate(double expected_rt, target_peptide, double maximal_rt_diff,
        libcpp_map[libcpp_string, int] already_seen, mst_settings * settings):
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
    # Select matching peakgroups from the target run (within the user-defined maximal rt deviation)
    cdef libcpp_pair[double, c_peakgroup*] retval
    retval.second = NULL

    cdef libcpp_vector[c_peakgroup].iterator pg_it
    cdef libcpp_vector[c_peakgroup*] matching_pg 
    cdef CyPrecursorWrapperOnly cypr
    cdef c_precursor * pep
    cdef libcpp_string id_string

    for precursor in target_peptide.getAllPrecursors():
        cypr = precursor
        pg_it = cypr.get_all_peakgroups_cy_begin()
        pep = cypr.getPeptidePrecursor()
        while pg_it != cypr.get_all_peakgroups_cy_end():
            if abs(deref(pg_it).normalized_retentiontime - expected_rt) < maximal_rt_diff:
                if deref(pg_it).fdr_score < settings.aligned_fdr_cutoff:
                    id_string = deref(pg_it).internal_id_ + deref(pep).curr_id_
                    if already_seen.find(id_string) == already_seen.end():
                        matching_pg.push_back( address(deref(pg_it) ))
            inc(pg_it)

    # If there are no peak groups present in the target run, we simply
    # return the expected retention time.
    ### if len(matching_peakgroups) == 0:
    if matching_pg.empty():
        retval.first = expected_rt
        return retval
        # return None, expected_rt

    # Select best scoring peakgroup among those in the matching RT window
    ### bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))
    cdef libcpp_vector[ c_peakgroup* ].iterator iter_m = matching_pg.begin()
    cdef libcpp_vector[ c_peakgroup* ].iterator bestScoring_pg = matching_pg.begin()
    cdef double best_fdr = deref(iter_m).fdr_score
    cdef int nr_below_fdr = 0
    cdef int nr_below_align = 0
    while iter_m != matching_pg.end():
        if (deref(iter_m).fdr_score <= best_fdr):
            best_fdr = deref(iter_m).fdr_score
            bestScoring_pg = iter_m
        if (deref(iter_m).fdr_score <= settings.fdr_cutoff):
            nr_below_fdr += 1
        if (deref(iter_m).fdr_score <= settings.aligned_fdr_cutoff):
            nr_below_align += 1
        inc(iter_m)

    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._aligned_fdr_cutoff]) > 1:
    ###     self.nr_multiple_align += 1
    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._fdr_cutoff]) > 1:
    ###     self.nr_ambiguous += 1
    if (nr_below_fdr > 1):
            settings.nr_ambiguous += 1
    if (nr_below_align > 1):
            settings.nr_multiple_align += 1

    # Decide which retention time to return:
    #  - the threading one based on the alignment
    #  - the one of the best peakgroup
    retval.second = deref(bestScoring_pg)
    if settings.correctRT_using_pg:
        retval.first = deref(bestScoring_pg).normalized_retentiontime
        return retval
    else:
        retval.first = expected_rt
        return retval

