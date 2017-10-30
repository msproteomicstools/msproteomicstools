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

include "PeakgroupWrapper.pyx"
include "PrecursorWrapper.pyx"
include "PrecursorGroup.pyx"

cdef class CyPrecursor(object):
    """ A set of peakgroups that belong to the same precursor in a single run.

    Each precursor has a backreference to its precursor group (heavy/light
    pair) it belongs to, the run it belongs to as well as its amino acid sequence.
    Furthermore, a unique id for the precursor and the protein name are stored.

    A precursor can return its best transition group, the selected peakgroup,
    or can return the transition group that is closest to a given iRT time.
    Its id is the transition_group_id (e.g. the id of the chromatogram)

    The "selected" peakgroup is represented by the peakgroup that belongs to
    cluster number 1 (cluster_id == 1) which in this case is "special".

    == Implementation details ==
    
    For memory reasons, we store all information about the peakgroup in a
    tuple (invariable). This tuple contains a unique feature id, a score and
    a retention time. Additionally, we also store, in which cluster the
    peakgroup belongs (if the user sets this).

    A peakgroup has the following attributes: 
        - an identifier that is unique among all other precursors 
        - a set of peakgroups 
        - a back-reference to the run it belongs to
    """

    cdef bool _decoy
    cdef libcpp_vector[c_peakgroup] cpeakgroups_ 
    cdef libcpp_string curr_id_
    cdef libcpp_string protein_name_ 
    cdef libcpp_string sequence_ 
    cdef CyPrecursorGroup precursor_group
    cdef object run

    def __init__(self, bytes this_id, run):
        self.curr_id_ = libcpp_string(<char*> this_id)
        self.run = run
        self._decoy = False

        # These remain NULL / unset:
        # cdef libcpp_vector[c_peakgroup] cpeakgroups_ 
        # cdef libcpp_string protein_name_ 
        # cdef libcpp_string sequence_ 
        # cdef object precursor_group
  
    def set_precursor_group(self, object p):
        self.precursor_group = p

    def set_decoy(self, bytes decoy):
        if decoy in ["FALSE", "False", "0"]:
            self._decoy = False
        elif decoy in ["TRUE", "True", "1"]:
            self._decoy = True
        else:
            raise Exception("Unknown decoy classifier '%s', please check your input data!" % decoy)

    def get_decoy(self):
        return self._decoy

    def setSequence(self, bytes p):
        self.sequence_ = libcpp_string(<char*> p)

    def getSequence(self):
        return <bytes>( self.sequence_ )

    def getProteinName(self):
        return <bytes>( self.protein_name_ )

    def setProteinName(self, bytes p):
        self.protein_name_ = libcpp_string(<char*> p)

    def get_id(self):
        return <bytes>( self.curr_id_ )
  
    def getPrecursorGroup(self):
        raise Exception("Wont do that! ")
        return self.precursor_group 

    def getPrecursorGroupId(self):
        raise Exception("Wont do that! ")
        return self.precursor_group.getPeptideGroupLabel()

    cdef bytes getRunId(self):
        return <bytes>self.run.get_id()

    def getRunId(self):
        return self.run.get_id()

    def getRun(self):
        raise Exception("Not doing that")
        return self.run

    def get_run_id(self):
      raise Exception("Not implemented")

    def __str__(self):
        return "%s (run %s)" % (self.get_id(), self.run)

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
        self.cpeakgroups_.push_back(pg)

    # 
    # Peakgroup cluster membership
    # 
    def select_pg(self, bytes this_id):
        self._setClusterID(this_id, 1)

    def unselect_pg(self, bytes this_id):
        self._setClusterID(this_id, -1)

    def setClusterID(self, bytes this_id, int cl_id):
        self._setClusterID(this_id, cl_id)

    cdef _setClusterID(self, bytes this_id, int cl_id):
        cdef libcpp_string s
        s = libcpp_string(<char*> this_id)
        cdef libcpp_vector[c_peakgroup].iterator it = self.cpeakgroups_.begin()
        cdef int nr_hit = 0
        while it != self.cpeakgroups_.end():
            if deref(it).internal_id_ == s:
                deref(it).cluster_id_ = cl_id
                nr_hit += 1
            inc(it)

        if nr_hit != 1:
              raise Exception("Error, found more than one peakgroup")

    def unselect_all(self):
        cdef libcpp_vector[c_peakgroup].iterator it = self.cpeakgroups_.begin()
        while it != self.cpeakgroups_.end():
            deref(it).cluster_id_ = -1
            inc(it)

    # 
    # Peakgroup selection
    # 
    def get_best_peakgroup(self):
        """
        Python code:
        ### if len(self.peakgroups_) == 0:
        ###     return None

        ### best_score = self.peakgroups_[0][1]
        ### result = self.peakgroups_[0]
        ### for peakgroup in self.peakgroups_:
        ###     if peakgroup[1] <= best_score:
        ###         best_score = peakgroup[1]
        ###         result = peakgroup
        ### index = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == result[0]][0]
        ### return MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])
        """
        if self.cpeakgroups_.empty():
             return None

        cdef libcpp_vector[c_peakgroup].iterator it = self.cpeakgroups_.begin()
        cdef libcpp_vector[c_peakgroup].iterator best = self.cpeakgroups_.begin()
        cdef double best_score = deref(it).fdr_score
        while it != self.cpeakgroups_.end():
            if deref(it).fdr_score <= best_score:
                best_score = deref(it).fdr_score
                best = it
            inc(it)

        result = CyPeakgroupWrapperOnly()
        result.inst = address(deref(best))
        result.peptide = self
        return result

    ### def _fixSelectedPGError(self, fixMethod="Exception"):
    ###   selected = [i for i,pg in enumerate(self.cluster_ids_) if pg == 1]
    ###   if len(selected) > 1:
    ###       print("Potential error detected in %s:\nWe have more than one selected peakgroup found. Starting error handling by using method '%s'." % (self, fixMethod))
    ###       best_score = self.peakgroups_[0][1]
    ###       best_pg = 0
    ###       for s in selected:
    ###           if best_score > self.peakgroups_[s][1]:
    ###               best_score = self.peakgroups_[s][1]
    ###               best_pg = s

    ###       if fixMethod == "Exception":
    ###           raise Exception("More than one selected peakgroup found in %s " % self )
    ###       elif fixMethod == "BestScore":
    ###           # Deselect all, then select the one with the best score...
    ###           for s in selected:
    ###               self.cluster_ids_[s] = -1
    ###           self.cluster_ids_[best_pg] = 1


    def get_selected_peakgroup(self):
        """
          return the selected peakgroup of this precursor, we can only select 1 or
          zero groups per chromatogram!

        Python code:

          #### # return the selected peakgroup of this precursor, we can only select 1 or
          #### # zero groups per chromatogram!
          #### selected = [i for i,pg in enumerate(self.cluster_ids_) if pg == 1]
          #### assert len(selected) < 2
          #### if len(selected) == 1:
          ####   index = selected[0]
          ####   result = self.peakgroups_[index]
          ####   return MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])
          #### else: 
          ####     return None

        """
        if self.cpeakgroups_.empty():
             return None

        cdef libcpp_vector[c_peakgroup].iterator it = self.cpeakgroups_.begin()
        cdef libcpp_vector[c_peakgroup].iterator best = self.cpeakgroups_.begin()
        cdef int nr_hit = 0
        while it != self.cpeakgroups_.end():
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

        Python code:
          ## selected = [i for i,pg in enumerate(self.cluster_ids_) if pg != -1]
          ## for index in selected:
          ##   result = self.peakgroups_[index]
          ##   yield MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])
        """

        cdef libcpp_vector[c_peakgroup].iterator it = self.cpeakgroups_.begin()
        while it != self.cpeakgroups_.end():
            if (deref(it).cluster_id_ == -1): 
                inc(it)
                continue
            result = CyPeakgroupWrapperOnly()
            result.inst = address(deref(it))
            result.peptide = self
            yield result
            inc(it)

    cdef libcpp_vector[c_peakgroup] get_all_peakgroups_cy(self):
        return self.cpeakgroups_

    cdef libcpp_vector[c_peakgroup].iterator get_all_peakgroups_cy_begin(self):
        return self.cpeakgroups_.begin()

    cdef libcpp_vector[c_peakgroup].iterator get_all_peakgroups_cy_end(self):
        return self.cpeakgroups_.end()

    def get_all_peakgroups(self):
        """
        Python code:
          ## for index, result in enumerate(self.peakgroups_):
          ##     yield MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])
        """

        cdef libcpp_vector[c_peakgroup].iterator it = self.cpeakgroups_.begin()
        while it != self.cpeakgroups_.end():
            result = CyPeakgroupWrapperOnly()
            result.inst = address(deref(it))
            result.peptide = self
            yield result
            inc(it)

    def getAllPeakgroups(self):
        return self.get_all_peakgroups()

    cdef libcpp_vector[c_peakgroup] getPeakGroupsVector(self):
        return self.cpeakgroups_
