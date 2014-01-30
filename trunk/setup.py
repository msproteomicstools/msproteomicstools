#!/usr/bin/env python

import os
from distutils.core import setup

import fnmatch
all_scripts = []
for root, dirnames, filenames in os.walk('analysis'):
  for filename in fnmatch.filter(filenames, '*.py'):
      all_scripts.append(os.path.join(root, filename))

setup(name='msproteomicstools',
      version='0.1.0',
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
          'openswathgui': 'gui/openswath',
      },
      scripts=all_scripts,
      description='Tools for MS-based proteomics',
      long_description='msproteomicstools - python module for MS-based proteomics',
      url='https://code.google.com/p/msproteomicstools',
      license='Modified BSD',
      platforms='any that supports python 2.7',
      )

