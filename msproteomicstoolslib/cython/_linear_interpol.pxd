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

cdef extern from "_linear_interpol.h":
    cdef cppclass c_linear_interpolate:
        c_linear_interpolate(libcpp_vector[double] & x, libcpp_vector[double] & y, double abs_err)
        double predict(double xnew)

cdef class CyLinearInterpolateWrapper(object):
    """
    Cython wrapper around c_linear_interpolate

    Another smoother that interpolates between the given data points. It is
    fast because its written in C++.

    This class expectes already smoothed x,y data (e.g. computed using a lowess
    or spline smoothing) but for applying the transformation, new x-values will
    be requested. The corresponding y values will be calculated by
    interpolation.

    The class provides the following methods:
        - `def __init__(self, x, y, double abs_err)`: initialize with two vectors, x and y
        - `def predict(self, list xnew)`: predict for Python
        - `cdef double predict_cy(self, double xnew)`: predict for Cython (low overhead)

    """
    cdef c_linear_interpolate * inst 

    cdef double predict_cy(self, double xnew)

