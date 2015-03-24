# Installation #

## Binary ##

For the TAPIR software, please see our website at http://proteomics.ethz.ch/tapir/ where binaries are available for Windows and Mac OS X. For a source installation of the package, see below.

## Dependencies ##

For the TAPIR graphical user interface, you will need the following dependencies:

  * Python >= 2.5
  * pymzML (from https://github.com/hroest/pymzML/archive/Release_0.7.5.zip)
  * numpy / scipy
  * PyQt4 >= 4.3
  * PyQwt >= 5.2
  * guidata
  * guiqwt

for linux, you can install those (except pymzml) on debian using `sudo apt-get install python python-numpy python-qt4 python-qwt5-qt4 python-guiqwt`

for Windows, you can download precompiled binaries, for example from http://www.lfd.uci.edu/~gohlke/pythonlibs/

pymzML needs to be installed separately, please download all necessary files (https://github.com/hroest/pymzML/archive/Release_0.7.5.zip), unzip and install with `python setup.py install`

## Install ##

After installing the dependencies, you can proceed to install msproteomicstools itself:

```
$ pip install msproteomicstools
```

Alternatively, a source install can be chosen as well:
```
$ svn checkout http://msproteomicstools.googlecode.com/svn/trunk/ msproteomicstools
$ python setup.py install --prefix=/your/install/path 
```

When performing a source install, the following packages need to be present (see setup.py):
  * numpy
  * scipy
  * cluster
  * pyteomics
  * xlsxwriter
  * xlwt
  * scikits.datasmooth
  * lxml



