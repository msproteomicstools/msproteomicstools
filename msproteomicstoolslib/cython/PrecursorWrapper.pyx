# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np

cdef extern from "precursor.h":
    cdef cppclass c_precursor:
        c_precursor()
        c_precursor(libcpp_string my_id, libcpp_string run_id)

        bool decoy
        libcpp_vector[c_peakgroup] peakgroups
        libcpp_string curr_id_
        libcpp_string protein_name_
        libcpp_string sequence_
        libcpp_string run_id_
        libcpp_string precursor_group_id

        libcpp_string getRunId()
        libcpp_string get_id()

        void add_peakgroup_tpl(c_peakgroup & pg, libcpp_string tpl_id, int cluster_id)


cdef class CyPrecursorWrapperOnly(object):
    """ A set of peakgroups that belong to the same precursor in a single run.

    Each precursor has a backreference to its precursor group identifier it
    belongs to, the run it belongs to as well as its amino acid sequence and
    protein name.

    Each precursor has a list of :class:`.CyPeakgroupWrapperOnly` that are
    found in the chromatogram of this precursor in this particular run.
    """

    cdef c_precursor * inst 
    cdef bool own_ptr

    def __dealloc__(self):
        if self.own_ptr:
            del self.inst

    def __init__(self, bytes this_id, run, take_ownership=True):
        if run is not None:
            self.inst = new c_precursor(this_id, <bytes>run.get_id())

        if take_ownership: self.own_ptr = True
        else: self.own_ptr = False

    #
    ### Getters
    cdef libcpp_string get_id_c(self):
        return deref(self.inst).curr_id_ 

    def get_id(self):
        return <bytes>( deref(self.inst).curr_id_ )

    def get_decoy(self):
        return ( deref(self.inst).decoy )

    def getSequence(self):
        return <bytes>( deref(self.inst).sequence_ )

    def getProteinName(self):
        return <bytes>( deref(self.inst).protein_name_ )

    def getRunId(self):
        return <bytes>( deref(self.inst).run_id_ )

    cdef libcpp_vector[c_peakgroup].iterator get_all_peakgroups_cy_begin(self):
        return deref(self.inst).peakgroups.begin()

    cdef libcpp_vector[c_peakgroup].iterator get_all_peakgroups_cy_end(self):
        return deref(self.inst).peakgroups.end()

    cdef c_precursor * getPeptidePrecursor(self):
        return self.inst

    #
    ### Setters
    def setProteinName(self, bytes p):
        deref(self.inst).protein_name_ = libcpp_string(<char*> p)

    def setSequence(self, bytes p):
        deref(self.inst).sequence_ = libcpp_string(<char*> p)

    def set_decoy(self, bytes decoy):
        if decoy in ["FALSE", "False", "0"]:
            deref(self.inst).decoy = False
        elif decoy in ["TRUE", "True", "1"]:
            deref(self.inst).decoy = True
        else:
            raise Exception("Unknown decoy classifier '%s', please check your input data!" % decoy)

    def set_precursor_group(self, precursor_group):
        deref(self.inst).precursor_group_id = <bytes>precursor_group.getPeptideGroupLabel()

    def add_peakgroup_tpl(self, pg_tuple, bytes tpl_id, int cluster_id=-1):
        """Adds a peakgroup to this precursor.

        The peakgroup should be a tuple of length 4 with the following components:
            0. id
            1. quality score (FDR)
            2. retention time (normalized)
            3. intensity
            (4. d_score optional)
        """
        # Check that the peak group is added to the correct precursor
        if self.get_id() != tpl_id:
            raise Exception("Cannot add a tuple to this precursor with a different id")

        if len(pg_tuple) == 4:
            pg_tuple = pg_tuple + (None,)

        assert len(pg_tuple) == 5
        cdef c_peakgroup pg
        pg.fdr_score = pg_tuple[1]
        pg.normalized_retentiontime = pg_tuple[2]
        pg.intensity_ = pg_tuple[3]
        pg.dscore_ = pg_tuple[4]
        pg.cluster_id_ = cluster_id
        pg.internal_id_ = libcpp_string(<char*> pg_tuple[0])
        pg.precursor = self.inst

        deref(self.inst).add_peakgroup_tpl(pg, tpl_id, cluster_id)


    # 
    # Peakgroup selection
    # 
    def get_best_peakgroup(self):
        """
        """
        if deref(self.inst).peakgroups.empty():
             return None

        cdef libcpp_vector[c_peakgroup].iterator it = deref(self.inst).peakgroups.begin()
        cdef libcpp_vector[c_peakgroup].iterator best = deref(self.inst).peakgroups.begin()
        cdef double best_score = deref(it).fdr_score
        while it != deref(self.inst).peakgroups.end():
            if deref(it).fdr_score <= best_score:
                best_score = deref(it).fdr_score
                best = it
            inc(it)

        result = CyPeakgroupWrapperOnly()
        result.inst = address(deref(best))
        result.peptide = self
        return result

    def unselect_all(self):
        cdef libcpp_vector[c_peakgroup].iterator it = deref(self.inst).peakgroups.begin()
        while it != deref(self.inst).peakgroups.end():
            deref(it).cluster_id_ = -1
            inc(it)

    def getAllPeakgroups(self):
        return self.get_all_peakgroups()

    def get_all_peakgroups(self):
        """
        """

        cdef libcpp_vector[c_peakgroup].iterator it = deref(self.inst).peakgroups.begin()
        while it != deref(self.inst).peakgroups.end():
            result = CyPeakgroupWrapperOnly()
            result.inst = address(deref(it))
            result.peptide = self
            yield result
            inc(it)

    def get_selected_peakgroup(self):
        """
          return the selected peakgroup of this precursor, we can only select 1 or
          zero groups per chromatogram!
        """
        ### print ("wrapper only - def get_selected_peakgroup(self):")
        ### if self.inst == NULL:
        ###     print ("is null!")
        ### else:
        ###     print ("is not null!")

        if deref(self.inst).peakgroups.empty():
             return None

        ## Somewhere here the segfault happens because peakgroups is not instantiated
        ## -> peakgroups.size() is some crazy number!
        ## print ("wrapper only - start: get_selected_peakgroup(self):")
        ## print ("size",  deref(self.inst).peakgroups.size())
        cdef libcpp_vector[c_peakgroup].iterator it = deref(self.inst).peakgroups.begin()
        cdef libcpp_vector[c_peakgroup].iterator best = deref(self.inst).peakgroups.begin()
        cdef int nr_hit = 0
        while it != deref(self.inst).peakgroups.end():
            if deref(it).cluster_id_ == 1:
                best = it
                nr_hit += 1
            inc(it)

        if nr_hit > 1:
              raise Exception("Error, found more than one peakgroup")
        if nr_hit == 0:
            return None

        result = CyPeakgroupWrapperOnly()
        result.inst = address(deref(best))
        result.peptide = self
        return result

    def getClusteredPeakgroups(self):
        """
        """

        cdef libcpp_vector[c_peakgroup].iterator it = deref(self.inst).peakgroups.begin()
        while it != deref(self.inst).peakgroups.end():
            if (deref(it).cluster_id_ == -1): 
                inc(it)
                continue
            result = CyPeakgroupWrapperOnly()
            result.inst = address(deref(it))
            result.peptide = self
            yield result
            inc(it)

    def getAllPrecursors(self):
        # dummy method!
        return [self]

    def printAddresses(self):
        print ("memaddr: {0:x}".format(<long> self.inst) )
        cdef libcpp_vector[c_peakgroup].iterator it = deref(self.inst).peakgroups.begin()
        while it != deref(self.inst).peakgroups.end():
            print (deref(it).internal_id_, " cl ", deref(it).cluster_id_, "memaddr: {0:x}".format(<long> address(deref(it))) )
            inc(it)

