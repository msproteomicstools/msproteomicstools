# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np

from libcpp.string cimport string as libcpp_string
from libcpp.vector  cimport vector as libcpp_vector
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from libcpp cimport bool

"""
cython -a --cplus msproteomicstoolslib/_optimized.pyx &&  python setup.py  build && cp build/lib.linux-x86_64-2.7/msproteomicstoolslib/_optimized.so msproteomicstoolslib/
"""

import numpy as np
import operator

cdef class CyLightTransformationData(object):
    """
    A lightweight data structure to store a transformation between retention times of multiple runs.
    """

    cdef dict data 
    cdef dict trafo 
    cdef dict stdevs
    cdef object reference

    def __init__(self, ref=None):
        self.data = {} 
        self.trafo = {} 
        self.stdevs = {} 
        self.reference = ref

    def addTrafo(self, run1, run2, CyLinearInterpolateWrapper trafo, stdev=None):
      """
      Add transformation between two runs
      """
      d = self.trafo.get(run1, {})
      d[run2] = trafo
      self.trafo[run1] = d

      d = self.stdevs.get(run1, {})
      d[run2] = stdev
      self.stdevs[run1] = d

    def addData(self, run1, data1, run2, data2, doSort=True):
      """
      Add raw data for the transformation between two runs
      """
      # Add data from run1 -> run2 and also run2 -> run1
      assert len(data1) == len(data2)
      self._doAddData(run1, data1, run2, data2, doSort)
      self._doAddData(run2, data2, run1, data1, doSort)

    def _doAddData(self, run1, data1, run2, data2, doSort):
      if doSort and len(data1) > 0:
          data1, data2 = zip(*sorted(zip(data1, data2)))
      data1 = np.array(data1)
      data2 = np.array(data2)
      d = self.data.get(run1, {})
      d[run2] = (data1,data2)
      self.data[run1] = d

    #
    ## Getters
    #
    def getData(self, run1, run2):
        return self.data[run1][run2]

    cdef CyLinearInterpolateWrapper getTrafoCy(self, run1, run2):
        return self.trafo[run1][run2]

    def getTrafo(self, run1, run2):
        return self.trafo[run1][run2]

    def getStdev(self, run1, run2):
        return self.stdevs[run1][run2]

    cdef double getStdevCy(self, run1, run2):
        return <double>self.stdevs[run1][run2]

    def getTransformation(self, run1, run2):
        return self.trafo[run1][run2]

    def getReferenceRunID(self):
        return self.reference

ctypedef np.float32_t DATA_TYPE

cdef extern from "peakgroup.h":
    cdef cppclass c_peakgroup:
        c_peakgroup()

        double fdr_score 
        double normalized_retentiontime 
        libcpp_string internal_id_ 
        double intensity_
        double dscore_ 
        int cluster_id_

cdef extern from "peakgroup.h":
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

        void add_peakgroup_tpl(c_peakgroup & pg, libcpp_string tpl_id, int cluster_id)

cdef extern from "peakgroup.h":
    cdef cppclass c_linear_interpolate:
        c_linear_interpolate(libcpp_vector[double] & x, libcpp_vector[double] & y, double abs_err)
        double predict(double xnew)

cdef class CyLinearInterpolateWrapper(object):
    cdef c_linear_interpolate * inst 

    def __dealloc__(self):
        del self.inst

    def __init__(self, x, y, double abs_err):
        cdef libcpp_vector[double] v1 = x
        cdef libcpp_vector[double] v2 = y

        self.inst = new c_linear_interpolate(v1, v2, abs_err)

    def predict(self, list xnew):
        ynew = []
        cdef double x
        for xn in xnew:
            ynew.append( <double>deref(self.inst).predict( <double>xn) )
        return ynew

    cdef double predict_cy(self, double xnew):
        return deref(self.inst).predict(xnew)

cdef class CyPeakgroupWrapperOnly(object):
    """
    """

    cdef c_peakgroup * inst 
    # cdef CyPrecursor peptide
    cdef CyPrecursorWrapperOnly peptide

    def __dealloc__(self):
        pass

    def __init__(self):
        pass

    # Do not allow setting of any parameters (since data is not stored here)
    def set_fdr_score(self, fdr_score):
        raise Exception("Cannot set in immutable object")

    def set_normalized_retentiontime(self, normalized_retentiontime):
        raise Exception("Cannot set in immutable object")

    def set_feature_id(self, id_):
        raise Exception("Cannot set in immutable object")

    def set_intensity(self, intensity):
        raise Exception("Cannot set in immutable object")

    cdef CyPrecursor getPeptide(self):
        return self.peptide

    def getPeptide(self):
        return self.peptide

    def get_dscore(self):
        return self.dscore_

    ## Select / De-select peakgroup
    def select_this_peakgroup(self):
        ## self.peptide.select_pg(self.get_feature_id())
        deref(self.inst).cluster_id_ = 1

    ## Select / De-select peakgroup
    def setClusterID(self, id_):
        raise Exception("Not implemented!")
        ### self.cluster_id_ = id_
        ### self.peptide.setClusterID(self.get_feature_id(), id_)

    def get_cluster_id(self):
        return deref(self.inst).cluster_id_
  
    def __str__(self):
        return "PeakGroup %s at %s s in %s with score %s (cluster %s)" % (self.get_feature_id(), self.get_normalized_retentiontime(), "run?", self.get_fdr_score(), self.get_cluster_id())

    def get_value(self, value):
        raise Exception("Needs implementation")

    def set_value(self, key, value):
        raise Exception("Needs implementation")

    def set_fdr_score(self, fdr_score):
        self.fdr_score = fdr_score

    def get_fdr_score(self):
        return deref(self.inst).fdr_score

    def set_normalized_retentiontime(self, float normalized_retentiontime):
        self.normalized_retentiontime = normalized_retentiontime

    def get_normalized_retentiontime(self):
        return deref(self.inst).normalized_retentiontime

    cdef get_normalized_retentiontime_cy(self):
        return deref(self.inst).normalized_retentiontime

    def set_feature_id(self, id_):
        raise Exception("Not implemented!")
        # self.id_ = id_

    def get_feature_id(self):
        return <bytes>(deref(self.inst).internal_id_)

cdef class CyPrecursorWrapperOnly(object):

    cdef c_precursor * inst 

    def __dealloc__(self):
        del self.inst

    def __init__(self, bytes this_id, run):
        self.inst = new c_precursor(this_id, <bytes>run.get_id())

    #
    ### Getters
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
        if deref(self.inst).peakgroups.empty():
             return None

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
        return self.precursor_group 

    def getPrecursorGroupId(self):
        return self.precursor_group.getPeptideGroupLabel()

    cdef bytes getRunId(self):
        return <bytes>self.run.get_id()

    def getRunId(self):
        return self.run.get_id()

    def getRun(self):
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
  
cdef class CyPrecursorGroup(object):
    """A set of precursors that are isotopically modified versions of each other.

    A collection of precursors that are isotopically modified versions of the
    same underlying peptide sequence. Generally these are heavy/light forms.
    """

    __slots__ = ["peptide_group_label_", "run_", "precursors_"]

    cdef libcpp_string peptide_group_label_ 
    cdef object run_
    cdef list precursors_ # list of peptide precursors (CyPrecursor)

    def __init__(self, peptide_group_label, run):
        self.peptide_group_label_ = peptide_group_label  
        self.run_ = run
        self.precursors_ = []

    def __str__(self):
        return "PrecursorGroup %s" % (self.getPeptideGroupLabel())

    ### def __lt__(self, other):

    ###     if self.run_.get_id() == other.run_.get_id():
    ###         return self.getPeptideGroupLabel() > other.getPeptideGroupLabel()
    ###     else:
    ###         return self.run_.get_id() > other.run_.get_id()

    def __iter__(self):
        for precursor in self.precursors_:
            yield precursor

    def __classInvariant__(self):
        if len(self.precursors_) > 0:
            # All precursor sequences should all be equal to the first sequence
            assert(all( [precursor.getSequence() == self.precursors_[0].getSequence() for precursor in self.precursors_] )) 
        return True

    # @class_invariant(__classInvariant__)
    def getPeptideGroupLabel(self):
        """
        getPeptideGroupLabel(self)
        Get peptide group label
        """
        return self.peptide_group_label_
  
    # @class_invariant(__classInvariant__)
    def addPrecursor(self, CyPrecursor precursor):
        """
        addPrecursor(self, precursor)
        Add precursor to peptide group
        """
        precursor.set_precursor_group( self )
        self.precursors_.append(precursor)

    # @class_invariant(__classInvariant__)
    def getPrecursor(self, curr_id):
        """
        getPrecursor(self, curr_id)
        Get the precursor for the given transition group id
        """
        for precursor in self:
            if precursor.get_id() == curr_id:
                return precursor
        return None

    # @class_invariant(__classInvariant__)
    def getAllPrecursors(self):
        """
        getAllPrecursors(self)
        Return a list of all precursors in this precursor group
        """
        return list(self)

    # @class_invariant(__classInvariant__)
    def getAllPeakgroups(self):
        """
        getAllPeakgroups(self)
        Generator of all peakgroups attached to the precursors in this group
        """
        for pr in self.precursors_:
            for pg in pr.get_all_peakgroups():
                yield pg

    # @class_invariant(__classInvariant__)
    def getOverallBestPeakgroup(self):
        """
        getOverallBestPeakgroup(self)
        Get the best peakgroup (by fdr score) of all precursors contained in this precursor group
        """
        allpg = list(self.getAllPeakgroups())
        if len(allpg) == 0:
            return None

        minscore = min([pg.get_fdr_score() for pg in allpg])
        return [pg for pg in allpg if pg.get_fdr_score() <= minscore][0]

    def get_decoy(self):
        """
        Whether the current peptide is a decoy or not

        Returns:
            decoy(bool): Whether the peptide is decoy or not
        """
        if len(self.precursors_) == 0:
            return False

        return self.precursors_[0].get_decoy()

@cython.boundscheck(False)
@cython.wraparound(False)
def static_cy_findBestPGFromTemplate(double expected_rt, CyPrecursorGroup target_peptide, double max_rt_diff,
        already_seen, double aligned_fdr_cutoff, double fdr_cutoff, correctRT_using_pg,
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
    matching_peakgroups = []
    cdef CyPeakgroupWrapperOnly pg_wr
    cdef libcpp_vector[ c_peakgroup ] matching_v
    for pg_ in target_peptide.getAllPeakgroups():
        pg_wr = <CyPeakgroupWrapperOnly> pg_
        if (abs(float(pg_wr.get_normalized_retentiontime()) - float(expected_rt)) < max_rt_diff):
            if pg_wr.get_fdr_score() < aligned_fdr_cutoff:
                if pg_wr.get_feature_id() + pg_wr.getPeptide().get_id() not in already_seen:
                    matching_peakgroups.append(pg_wr)

    for pg_ in target_peptide.getAllPeakgroups():
        pg_wr = <CyPeakgroupWrapperOnly> pg_
        if (abs(float(deref(pg_wr.inst).normalized_retentiontime) - float(expected_rt)) < max_rt_diff):
            if deref(pg_wr.inst).fdr_score < aligned_fdr_cutoff:
                if deref(pg_wr.inst).internal_id_ + pg_wr.getPeptide().get_id() not in already_seen:
                    # matching_peakgroups.append(pg_wr)
                    matching_v.push_back( deref(pg_wr.inst) )

    #### matching_peakgroups = [pg_ for pg_ in target_peptide.getAllPeakgroups() 
    ####     if (abs(float(pg_.get_normalized_retentiontime()) - float(expected_rt)) < max_rt_diff) and
    ####         pg_.get_fdr_score() < aligned_fdr_cutoff and 
    ####         pg_.get_feature_id() + pg_.getPeptide().get_id() not in already_seen]

    #### cdef double rt
    #### cdef double fdr_score
    #### # matching_peakgroups = []
    #### cdef libcpp_vector[c_peakgroup] CPeakgroups =  target_peptide.getPeakGroupsVector()
    #### cdef libcpp_vector[c_peakgroup].iterator it = target_peptide.getPeakGroupsVector().begin()

    #### cdef libcpp_vector[c_peakgroup] matching_peakgroups

    #### while it != target_peptide.getPeakGroupsVector().end():
    #### # for pg_ in target_peptide.getAllPeakgroups():
    ####     ## pg_ = deref(it)
    ####     rt = deref(it).normalized_retentiontime
    ####     fdr_score = deref(it).fdr_score
    ####     if fdr_score < aligned_fdr_cutoff:
    ####         if (abs(rt - float(expected_rt)) < max_rt_diff):
    ####             if deref(it).internal_id_ + target_peptide.get_id() not in already_seen:
    ####                 matching_peakgroups.push_back( deref(it) )
    ####     inc(it)


    # If there are no peak groups present in the target run, we simply
    # return the expected retention time.
    ### if len(matching_peakgroups) == 0:
    if matching_v.empty():
        return None, expected_rt

    cdef libcpp_vector[ c_peakgroup ].iterator iter_m = matching_v.begin()
    cdef libcpp_vector[ c_peakgroup ].iterator bestScoring_pg = matching_v.begin()
    cdef double best_fdr = deref(iter_m).fdr_score
    while iter_m != matching_v.end():
        if (deref(iter_m).fdr_score <= best_fdr):
            best_fdr = deref(iter_m).fdr_score
            bestScoring_pg = iter_m
        inc(iter_m)


    # Select best scoring peakgroup among those in the matching RT window
    ### bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))

    bestScoringPG = CyPeakgroupWrapperOnly()
    bestScoringPG.inst = address(deref(bestScoring_pg))
    # bestScoringPG.peptide = target_peptide ## TODO:

    #### # Printing for debug mode
    #### if verbose and False:
    ####     closestPG = min(matching_peakgroups, key=lambda x: abs(float(x.get_normalized_retentiontime()) - expected_rt))
    ####     print("    closest:", closestPG.print_out(), "diff", abs(closestPG.get_normalized_retentiontime() - expected_rt) )
    ####     print("    bestScoring:", bestScoringPG.print_out(), "diff", abs(bestScoringPG.get_normalized_retentiontime() - expected_rt) )
    ####     print()

    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._aligned_fdr_cutoff]) > 1:
    ###     self.nr_multiple_align += 1
    ### if len([pg_ for pg_ in matching_peakgroups if pg_.get_fdr_score() < self._fdr_cutoff]) > 1:
    ###     self.nr_ambiguous += 1

    # Decide which retention time to return:
    #  - the threading one based on the alignment
    #  - the one of the best peakgroup
    if correctRT_using_pg:
        return bestScoringPG, deref(bestScoring_pg).normalized_retentiontime
    else:
        return bestScoringPG, expected_rt

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
cdef static_findBestPGFromTemplate(double expected_rt, CyPrecursorGroup target_peptide, double max_rt_diff,
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



