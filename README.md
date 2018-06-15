[![Build Status](https://travis-ci.org/msproteomicstools/msproteomicstools.svg?branch=master)](https://travis-ci.org/msproteomicstools/msproteomicstools) [![Project Stats](https://www.openhub.net/p/msproteomicstools/widgets/project_thin_badge.gif)](https://www.openhub.net/p/msproteomicstools)

This repository contains some hopefully useful tools for mass spectrometry
applied to proteomics. Speficically, it contains 

- the TAPIR visualization software (available in ./gui)
- the MS Proteomics Tools Library (available in ./msproteomicstoolslib)
- a set of executables and scripts (available under ./analysis) including the TRIC alignment tool

The code is under the 3-clause BSD licence (see the [LICENSE](LICENSE)
and the [AUTHORS](AUTHORS.txt)  files). For full documentation, see [the online
documentation](http://msproteomicstools.hroest.ch/index.html).

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

# MS Proteomics Tools Library 

The mass spectrometric (MS) Proteomics Tools Library contains multiple Python
functions useful for MS-based proteomics.

# Install

After installing the dependencies, you can proceed to install msproteomicstools itself:

    pip install numpy
    pip install msproteomicstools

If you are on windows, there is a chance the above will not work as Biopython
and numpy need to be compiled from source. Currently, the way to get these is
by installing [Anaconda](https://www.continuum.io/downloads) and then running
in the Anaconda prompt:

    conda install biopython
    pip install msproteomicstools

Alternatively, you can also install from source:

    git checkout https://github.com/msproteomicstools/msproteomicstools.git
    pip install numpy
    python setup.py install --prefix=/your/install/path 


## Optional packages 

There are some extra packages that can increase the features and improve the speed of the toolset

### Fast lowess

To get fast lowess performance (several orders of magnitude faster), do the
following

    git clone https://github.com/carljv/Will_it_Python.git
    cd Will_it_Python/MLFH/CH2/lowess\ work/
    python setup.py build
    sudo python setup.py install

### Rpy2

If you would like to use the rpy2 bridge, you should install the `rpy2` package.

## Documentation

Documentation for the library can be found at
[http://msproteomicstools.roestlab.org/](http://msproteomicstools.roestlab.org/) 
which contains source code documentation of the available functions and
objects.

