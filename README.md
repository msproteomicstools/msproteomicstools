This repository contains some hopefully useful tools for mass spectrometry
applied to proteomics. Speficically, it contains 

- the TAPIR visualization software (available in ./gui)
- the MS Proteomics Tools Library (available in ./msproteomicstoolslib)
- a set of executables and scripts (available under ./analysis)

The code is under the 3-clause BSD licence (see the [COPYRIGHT](COPYRIGHT.txt)
and the [AUTHORS](AUTHORS.txt)  files).

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


## Extra packages 

There are some extra packages that can increase the features and improve the speed of the toolset

### Fast lowess

To get fast lowess performance (several orders of magnitude faster), do the
following

    git clone https://github.com/carljv/Will_it_Python.git
    cd Will_it_Python/MLFH/CH2/lowess\ work/
    python setup.py build
    sudo python setup.py install


