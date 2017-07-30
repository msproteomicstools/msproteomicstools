# distutils: language = c++
# cython: c_string_type=str, c_string_encoding=ascii
cimport cython
cimport libc.stdlib
cimport numpy as np

from libcpp.string cimport string as libcpp_string

cdef class CyPrecursorGroup(object):
    """See :class:`.PrecursorGroup` for a description.

    This implementation is pure Cython.

    Attributes:
        - self.peptide_group_label_: Identifier or precursor group 
        - self.run_: Reference to the :class:`.Run` where this PrecursorGroup is from
        - self.precursors_: List of :class:`.CyPrecursorWrapperOnly`
    """

    cdef libcpp_string peptide_group_label_ 
    cdef object run_
    cdef list precursors_ # list of peptide precursors (CyPrecursor)

