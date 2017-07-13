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

    cdef CyLinearInterpolateWrapper getTrafoCy(self, libcpp_string run1, libcpp_string run2):
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

