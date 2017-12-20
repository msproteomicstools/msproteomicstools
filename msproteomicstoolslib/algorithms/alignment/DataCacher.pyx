# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np
import numpy as np

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
    """
    Wrapper around c_data_cache which allows storage of multiple lists rt/fdr values for later alignment.

    Basically the data cacher allows storage of an fdr and a rt vector for each of N runs.

    Retrieval occurs by asking for these vectors for a specific combination of
    runs, the cacher returns a vector of RT pairs that occur in both runs.
    """
    cdef c_data_cache * inst 

    def __dealloc__(self):
        del self.inst

    def __init__(self):
        self.inst = new c_data_cache()

    def appendValuesForPeptide(self, cached_values):
        """
        Store the retention time values for a single peptide across all runs.

        Parameters
        ---------
        cached_values : list( list( double ) )
            A list of length N (number of runs) where for each run a pair of
            values is provided: (fdr,rt). In case the peptide was identified,
            None can be provided instead.
        """
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

    def retrieveValues(self, int run1, int run2):
        """
        Retrieve all paired RT values for two given runs

        Parameters
        ---------
        run1 : int
            Index of the the first run
        run2 : int
            Index of the the second run


        Returns
        -------
            tuple(list(double), list(double)) : the two lists containing matched RT values

        """
        cdef libcpp_vector[double] rt1
        cdef libcpp_vector[double] rt2
        deref(self.inst).retrieveValues(rt1, rt2, run1, run2)
        ## takes a bit longer:
        ## return rt1, rt2

        cdef np.ndarray[np.float64_t, ndim=1] data1
        cdef np.ndarray[np.float64_t, ndim=1] data2

        cdef libcpp_vector[double].iterator it = rt1.begin()
        cdef int i = 0

        cdef unsigned int rt1_size = rt1.size()
        cdef unsigned int rt2_size = rt2.size()
        data1 = np.zeros( (rt1_size,), dtype=np.float64)
        data2 = np.zeros( (rt1_size,), dtype=np.float64)

        while it != rt1.end():
            data1[i] = deref(it)
            inc(it)
            i += 1

        it = rt2.begin()
        i = 0
        while it != rt2.end():
            data2[i] = deref(it)
            inc(it)
            i += 1

        return data1, data2

