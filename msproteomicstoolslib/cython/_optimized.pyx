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

###########################################################################
# Collection C++ includes
#  * peakgroup.h
#      Peakgroup model implemented in C++
#  * precursor.h
#      Precursor model implemented in C++
#  * _linear_interpol.h
#      Linear interpolation implemented in C++
#    
# Collection of Cython includes
# 
#  * _linear_interpol.pyx
#      Linear wrapper
#  * PeakgroupWrapper.pyx
#      Peakgroup wrapper
#  * PrecursorWrapper.pyx
#      Precursor wrapper
#  * LightTransformationData.pyx
#      Data structure for storing RT transformations
#  * PrecursorGroup.pyx
#      Data structure for storing precursor groups
#  
#  * MSTAlignment.pyx
#      Set of functions ported from the MST alignment code (in Python)
#  * MSTAlignment_fast.pyx
#      Set of functions ported from the MST alignment code (in Cython)
#

# note that Precursor.pyx is a free-standing module


include "LightTransformationData.pyx"
include "PeakgroupWrapper.pyx"
include "PrecursorWrapper.pyx"
include "PrecursorGroup.pyx"

include "MSTAlignment.pyx"
include "MSTAlignment_fast.pyx"

