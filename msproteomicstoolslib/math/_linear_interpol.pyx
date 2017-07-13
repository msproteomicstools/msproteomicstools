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

