# distutils: language = c++
# cython: c_string_type=str, c_string_encoding=ascii
# encoding: latin-1

import numpy as np
from _linear_interpol cimport CyLinearInterpolateWrapper, c_linear_interpolate
from libcpp.string cimport string as libcpp_string

cdef class CyLightTransformationData(object):
    """
    Cython implementation of :class:`.LightTransformationData`

    A lightweight data structure to store a transformation between retention times of multiple runs.
    """

    cdef dict data 
    cdef dict trafo 
    cdef dict stdevs
    cdef object reference

    cdef CyLinearInterpolateWrapper getTrafoCy(self, libcpp_string run1, libcpp_string run2)
    cdef double getStdevCy(self, run1, run2)

