#!/usr/bin/env python

import os
from distutils.core import setup

import fnmatch
all_scripts = []
for root, dirnames, filenames in os.walk('analysis'):
  for filename in fnmatch.filter(filenames, '*.py'):
      all_scripts.append(os.path.join(root, filename))
all_scripts.extend(["./gui/AlignmentGUI.py"])

setup(name='msproteomicstools',
      version='0.2.1',
      packages = ['msproteomicstoolslib', 
                  "msproteomicstoolslib.algorithms",
                  "msproteomicstoolslib.algorithms.alignment",
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
      scripts=all_scripts,
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
          "cluster >= 1.1.2",
          "pyteomics >= 2.4.0",
          "xlsxwriter >= 0.5.3 ",
          "numpy",
          "scipy",
          'xlwt',
          'scikits.datasmooth',
          'rpy2',
          'pymzml',
          'lxml'
      ],
      test_suite="nose.collector",
      tests_require="nose",
      )

