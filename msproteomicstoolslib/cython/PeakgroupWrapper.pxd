# distutils: language = c++
# cython: c_string_type=str, c_string_encoding=ascii
cimport cython
cimport libc.stdlib
cimport numpy as np

from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from libcpp cimport bool

from PrecursorWrapper cimport c_precursor, CyPrecursorWrapperOnly

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
    See :class:`.PeakGroupBase` for a detailed description.

    This implementation stores a pointer to a C++ object holding the actual
    data. The data access works very similarly as for any :class:`.PeakGroupBase`.
    """


    cdef c_peakgroup * inst 
    cdef CyPrecursorWrapperOnly peptide

    cdef get_normalized_retentiontime_cy(self)


