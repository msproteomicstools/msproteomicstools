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

cdef extern from "DataCacher.h":
    cdef cppclass c_data_cache:
        c_data_cache()
        void addValues(libcpp_vector[double] & fdr, libcpp_vector[double] & rt)
        void retrieveValues(libcpp_vector[double] & rt1, libcpp_vector[double] & rt2, int run1, int run2)

cdef class CyDataCacher(object):
    cdef c_data_cache * inst 

    def __dealloc__(self):
        del self.inst

    def __init__(self):
        self.inst = new c_data_cache()

    def appendValuesForPeptide(self, cached_values):
        cdef libcpp_vector[double] fdr
        cdef libcpp_vector[double] rt

        for cc in cached_values:
            if cc is None:
                fdr.push_back(-1)
                rt.push_back(-1)
            else:
                fdr_, rt_ = cc
                fdr.push_back(fdr_)
                rt.push_back(rt_)

        deref(self.inst).addValues(fdr, rt)

    def retrieveValues(self, data1, data2, int run1, int run2):
        cdef libcpp_vector[double] rt1
        cdef libcpp_vector[double] rt2
        deref(self.inst).retrieveValues(rt1, rt2, run1, run2)

        cdef libcpp_vector[double].iterator it = rt1.begin()
        while it != rt1.end():
            data1.append( deref(it) )
            inc(it)

        it = rt2.begin()
        while it != rt2.end():
            data2.append( deref(it) )
            inc(it)

