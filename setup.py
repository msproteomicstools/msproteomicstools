#!/usr/bin/env python

import os
import numpy
from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# get py-earth from https://github.com/jcrudy/py-earth/

import fnmatch
all_scripts = []
for root, dirnames, filenames in os.walk('analysis'):
  for filename in fnmatch.filter(filenames, '*.py'):
      all_scripts.append(os.path.join(root, filename))
all_scripts.extend(["./gui/AlignmentGUI.py"])
all_scripts.extend(["./gui/TAPIR.py"])

import sys
if (sys.version_info > (3, 0)):
    extra_installs = []
else:
    extra_installs = []


ext_modules = [
        cythonize(Extension('msproteomicstoolslib/cython/_optimized', sources=["msproteomicstoolslib/cython/_optimized.pyx"], language="c++", extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]))[0],
        cythonize("msproteomicstoolslib/cython/Precursor.pyx", language="c++")[0],
        cythonize("msproteomicstoolslib/cython/Precursor.pyx", language="c++")[0],
        cythonize("msproteomicstoolslib/algorithms/alignment/DataCacher.pyx", language="c++")[0]
        ]

setup(name='msproteomicstools',
      version='0.5.0',
      packages = ['msproteomicstoolslib', 
                  "msproteomicstoolslib.algorithms",
                  "msproteomicstoolslib.algorithms.alignment",
                  "msproteomicstoolslib.algorithms.shared",
                  "msproteomicstoolslib.algorithms.PADS",
                  "msproteomicstoolslib.data_structures",
                  "msproteomicstoolslib.format",
                  "msproteomicstoolslib.math",
                  "msproteomicstoolslib.util",
                  "openswathgui",
                  "openswathgui.models",
                  "openswathgui.views",
                 ],
      package_dir = {
          'openswathgui': 'gui/openswathgui',
      },
      package_data={'msproteomicstoolslib.data_structures':
          ['modifications_default.tsv']},
      scripts=all_scripts,
      ext_modules=ext_modules,
      include_dirs=[numpy.get_include()],
      description='Tools for MS-based proteomics',
      long_description='msproteomicstools - python module for MS-based proteomics',
      url='https://code.google.com/p/msproteomicstools',
      license='Modified BSD',
      platforms='any that supports python 2.7',
      classifiers=[
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: BSD License',
      'Operating System :: OS Independent',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Topic :: Scientific/Engineering :: Chemistry',
      ],
      install_requires=[
          "numpy",
          "scipy",
          "cluster == 1.2.2", # note that 1.1.2 does not work with py3
          "pyteomics >= 2.4.0",
          "statsmodels >= 0.6.0",
          "xlsxwriter >= 0.5.3 ", # for xlsx
          # 'xlwt', # for xls
          'scikits.datasmooth',
          # versions 7.6 and 7.7 are broken for us (use spectra sanity check)
          'pymzml == 0.7.5',
          'lxml',
          'configobj',
          'biopython',
          'xlwt',
      ] + extra_installs,
      extras_require = {
          'RSmoothing' : ["rpy2"]
      },
      test_suite="nose.collector",
      tests_require="nose",
      )

