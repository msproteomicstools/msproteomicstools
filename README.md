[![Build Status](https://travis-ci.org/msproteomicstools/msproteomicstools.svg?branch=master)](https://travis-ci.org/msproteomicstools/msproteomicstools) [![Project Stats](https://www.openhub.net/p/msproteomicstools/widgets/project_thin_badge.gif)](https://www.openhub.net/p/msproteomicstools)

This repository contains some hopefully useful tools for mass spectrometry
applied to proteomics. Speficically, it contains 

- the TAPIR visualization software (available in ./gui)
- the MS Proteomics Tools Library (available in ./msproteomicstoolslib)
- a set of executables and scripts (available under ./analysis) including the TRIC alignment tool

The code is under the 3-clause BSD licence (see the [LICENSE](LICENSE)
and the [AUTHORS](AUTHORS.txt)  files). For full documentation, see [the online
documentation](http://msproteomicstools.hroest.ch/index.html).

# Install

You can install the library and tools as follows:

    pip install numpy
    pip install pymzml==0.7.8
    pip install Biopython
    pip install msproteomicstools

### Fast lowess

To get fast lowess performance (several orders of magnitude faster), do the
following

    git clone https://github.com/carljv/Will_it_Python.git
    cd Will_it_Python/MLFH/CH2/lowess\ work/
    python setup.py build
    sudo python setup.py install

# TRIC

The ./analysis folder contains the TRIC alignment tool, see [the TRIC
manual](TRIC-README.md) for further information.

# TAPIR 

TAPIR is a software that allows visualization of targeted proteomics data. See
[the TAPIR install instructions](INSTALL-TAPIR.md) for further information on
installation.

Please also see the [project homepage](http://proteomics.ethz.ch/tapir/)
for further information.  Binaries are available for Windows and Mac OS X. For
a source installation of the package, see below.

# Executables

The ./analysis folder contains multiple potentially useful executables,
including a tools for high throughput targeted proteomics data analysis (such
as SWATH-MS data analysis). See [the TRIC manual](TRIC-README.md) for further
information on the TRIC alignment tool.

# Documentation

Documentation for the library can be found at
[http://msproteomicstools.roestlab.org/](http://msproteomicstools.roestlab.org) 
which contains source code documentation of the available functions and
objects.

