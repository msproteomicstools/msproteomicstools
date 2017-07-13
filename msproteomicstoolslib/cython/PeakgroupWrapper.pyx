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

cdef extern from "peakgroup.h":
    cdef cppclass c_peakgroup:
        c_peakgroup()

        double fdr_score 
        double normalized_retentiontime 
        libcpp_string internal_id_ 
        double intensity_
        double dscore_ 
        int cluster_id_
        c_precursor * precursor

        c_precursor * getPeptide()


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
