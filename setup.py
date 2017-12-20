#!/usr/bin/env python

import os
import numpy
from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# get py-earth from https://github.com/jcrudy/py-earth/

import sys
with_cython = False
if "--with_cython" in sys.argv:
    with_cython = True
    sys.argv.remove("--with_cython")

import fnmatch
all_scripts = []
for root, dirnames, filenames in os.walk('analysis'):
  for filename in fnmatch.filter(filenames, '*.py'):
      all_scripts.append(os.path.join(root, filename))
all_scripts.extend(["./gui/AlignmentGUI.py"])
all_scripts.extend(["./gui/TAPIR.py"])

extra_installs = []

ext_modules = []
if with_cython:
    ext_modules = [
	    cythonize(Extension('msproteomicstoolslib/cython/_optimized', sources=["msproteomicstoolslib/cython/_optimized.pyx"], language="c++", extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]))[0],
	    cythonize("msproteomicstoolslib/cython/Precursor.pyx", language="c++")[0],
	    cythonize("msproteomicstoolslib/algorithms/alignment/DataCacher.pyx", language="c++")[0]
	    ]

setup(name='msproteomicstools',
      version='0.7.0',
      description='Tools for MS-based proteomics',
      long_description='msproteomicstools - python module for MS-based proteomics',
      url='https://github.com/msproteomicstools/msproteomicstools',
      license='Modified BSD',
      platforms=["any"],
      classifiers=[
      'Environment :: Console',
      'Environment :: X11 Applications :: Qt',

      'Intended Audience :: Science/Research',
      'Intended Audience :: Developers',

      'License :: OSI Approved :: BSD License',

      'Topic :: Documentation :: Sphinx',

      'Operating System :: OS Independent',

      # Supported Python versions:
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python :: 3.4',
      'Programming Language :: Python :: 3.5',
      'Programming Language :: Python :: 3.6',

      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Topic :: Scientific/Engineering :: Chemistry',
      ],

      # Package and install info
      scripts=all_scripts,
      packages = ['msproteomicstoolslib', 
                  "msproteomicstoolslib.cython",
                  "msproteomicstoolslib.algorithms",
                  "msproteomicstoolslib.algorithms.alignment",
                  "msproteomicstoolslib.algorithms.shared",
                  "msproteomicstoolslib.algorithms.PADS",
                  "msproteomicstoolslib.algorithms.graphs",
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
      ext_modules=ext_modules,
      include_dirs=[numpy.get_include()],
      test_suite="nose.collector",
      tests_require="nose",
)

