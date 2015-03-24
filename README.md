This repository contains some hopefully useful tools for mass spectrometry
applied to proteomics. Speficically, it contains 

- the TAPIR visualization software
- the MS Proteomics Tools Library 
- a set of executables and scripts (available under ./analysis)

# TAPIR 

TAPIR is a software that allows visualization of targeted proteomics data. See
[the TAPIR install instructions](INSTALL-TAPIR) for further information on
installation.

Please also see the [project homepage](http://proteomics.ethz.ch/tapir/ TAPIR)
for further information.  Binaries are available for Windows and Mac OS X. For
a source installation of the package, see below.

# MS Proteomics Tools Library 

The mass spectrometric (MS) Proteomics Tools Library contains multiple Python
functions useful for MS-based proteomics.

## Documentation

Documentation for the library can be found at [http://proteomics.ethz.ch/msproteomicstools/ ](http://proteomics.ethz.ch/msproteomicstools/) which contains source code documentation of the available functions and objects.

## Install ##

After installing the dependencies, you can proceed to install msproteomicstools itself:

    pip install numpy
    pip install msproteomicstools

Alternatively, a source install can be chosen as well:

    git checkout https://github.com/msproteomicstools/msproteomicstools.git
    pip install numpy
    python setup.py install --prefix=/your/install/path 

## Testing 

After installation, you can run all tests 

    nosetests test

with coverage analysis

    nosetests test --with-coverage --cover-package=msproteomicstoolslib

or disable all slow tests for a faster response

    nosetests -a '!slow' test

## Dependencies ##

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


## Extra packages 

There are some extra packages that can increase the features and improve the speed of the toolset

### Fast lowess

To get fast lowess performance (several orders of magnitude faster), do the
following

git clone https://github.com/carljv/Will_it_Python.git
cd Will_it_Python/MLFH/CH2/lowess\ work/
python setup.py build
sudo python setup.py install


