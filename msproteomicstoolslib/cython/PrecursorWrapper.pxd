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

from PeakgroupWrapper cimport c_peakgroup

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

    #
    ### Getters
    cdef libcpp_string get_id_c(self)

    cdef libcpp_vector[c_peakgroup].iterator get_all_peakgroups_cy_begin(self)

    cdef libcpp_vector[c_peakgroup].iterator get_all_peakgroups_cy_end(self)

    cdef c_precursor * getPeptidePrecursor(self)

