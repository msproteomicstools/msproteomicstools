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
#  * MSTAlignment_fast.pyx
#      Set of functions ported from the MST alignment code (in Cython)
#

from PeakgroupWrapper cimport CyPeakgroupWrapperOnly
from PeakgroupWrapper cimport c_peakgroup

from PrecursorWrapper cimport CyPrecursorWrapperOnly
from PrecursorWrapper cimport c_precursor
from PrecursorGroup cimport CyPrecursorGroup

include "MSTAlignment_fast.pyx"

