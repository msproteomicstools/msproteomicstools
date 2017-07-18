Please follow the install instructions in the [main README file](README.md) before you continue reading

# Testing 

After installation, you can run all tests 

    nosetests test

with coverage analysis

    nosetests test --with-coverage --cover-package=msproteomicstoolslib

or disable all slow tests for a faster response

    nosetests -a '!slow' test

using the same approach, also tests that require rpy2 or cython can be disabled:

    nosetests -a '!slow,!rpy2,!cython' test

# Cython

There is an (optional) Cython extension that currently only works under Python 2.x, to activate it use

	python setup.py build --with_cython

the resulting C++ modules are intended to speed up TRIC substantially.

# Dependencies

For the TAPIR graphical user interface, you will need the following dependencies:

  * Python >= 2.5
  * pymzML >= 0.7.5
  * numpy / scipy
  * PyQt4 >= 4.3
  * PyQwt >= 5.2
  * guidata
  * guiqwt

To install on Linux, you can install those (except pymzml) on debian using:

    sudo apt-get install python python-numpy python-qt4 python-qwt5-qt4 python-guiqwt

For Windows, you can download precompiled binaries, for example from
http://www.lfd.uci.edu/~gohlke/pythonlibs/ to obtain the binary Python modules.

# Documentation 

To build the documentations, please see the [documentation README](docs/README)

# PyPI upload

To upload a new version to PyPI:

    python setup.py sdist register upload

