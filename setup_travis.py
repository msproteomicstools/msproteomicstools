#!/usr/bin/env python

import os
import numpy
from setuptools import setup
from Cython.Build import cythonize

import sys
with_cython = False
if "--with_cython" in sys.argv:
    with_cython = True
    sys.argv.remove("--with_cython")

# only need this to copy over the shared libraries
ext_modules = []
if with_cython:
    ext_modules = [
	    cythonize("msproteomicstoolslib/cython/_optimized.pyx", language="c++")[0],
	    cythonize("msproteomicstoolslib/cython/Precursor.pyx", language="c++")[0],
	    cythonize("msproteomicstoolslib/algorithms/alignment/DataCacher.pyx", language="c++")[0]
	    ]

setup(name='msproteomicstools',
      # Package and install info
      ext_modules=ext_modules,
      include_dirs=[numpy.get_include()]
)

