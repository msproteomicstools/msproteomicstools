
# General Information

This document contains information on how to install TAPIR from source. Binaries are available for Windows and Mac OSX at the
[project homepage](http://proteomics.ethz.ch/tapir/) which is what you should try first before proceeding here. 

Dependencies for TAPIR are the following:

- pymzML (from `https://github.com/hroest/pymzML/archive/Release_0.7.5.zip`)
- Python >= 2.5
- numpy / scipy
- PyQt4 >= 4.3
- PyQwt >= 5.2
- PIL (Python Imaging Library, maybe not necessary...)
- guidata
- guiqwt


# Linux install 

- Debian/Ubuntu/ArchLinux: guiqwt should be packaged, basically install these
  packages

```
sudo apt-get install python python-numpy \
     python-qt4 python-qwt5-qt4 python-guiqwt
```

- RedHat: PyQwt, PyQt, scipy, numpy, python-imaging should be the packages

Then install the pymzml and msproteomicstools:
```
git clone https://github.com/pymzml/pymzML.git
python setup.py install
svn checkout http://msproteomicstools.googlecode.com/svn/trunk/ msproteomicstools
cd msproteomicstool
```

# Windows

First install Python (e.g. version 2.7) from `http://www.python.org/getit/` or through Anaconda `https://www.continuum.io/downloads`

then go to the following website and download these binary packages (the
following assumes 32 bit window, for 64 bit please choose the appropriate
package):

`http://www.lfd.uci.edu/~gohlke/pythonlibs/`

1. Get numpy (numpy-MKL-1.7.2rc1.win32-py2.7.exe)
2. Get PyQt (PyQt-Py2.7-x32-gpl-4.9.6-1.exe)
3. Get PyQwt (PyQwt-5.2.1-py2.7-x32-pyqt4.9.6-numpy1.7.1.exe)
4. Get guiqwt (guiqwt-2.3.1.win32-py2.7.exe)
5. Get guidata (guidata-1.6.1.win32-py2.7.exe)

Next, open a command prompt (type `cmd.exe` in your search bar in Windows and click on it) and go to your home directory: 

`cd \my\home\dir`

If you do not have an SVN or git client, one essay solution on Windows is to
use `http://downloadsvn.codeplex.com/` to checkout the pymzml and
msproteomicstools projects.

Next, download and install pymzml using the following link:

    https://github.com/hroest/pymzML/archive/Release_0.7.5.zip

Then navigate to the download location:

`C:\Python27\python.exe setup.py install`

next, download msproteomicstools

```
svn checkout http://msproteomicstools.googlecode.com/svn/trunk/ msproteomicstools
cd msproteomicstool
```

then edit `gui\openswath\AlignmentGUI.py` and set USE_GUIQWT to False (replace
"True" with "False"). This should fix the script and you can start TAPIR

```C:\Python27\python.exe gui\openswath\AlignmentGUI.py```

# Mac OS X

Thanks to George Rosenberger

This install was tested on OS X 10.8

1. Install msproteomicstools

```bash
svn checkout http://msproteomicstools.googlecode.com/svn/trunk/ msproteomicstools-read-only
```


2. Install Homebrew

```bash
ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go)"
brew doctor
```

3. Dependencies

```bash
brew install python
brew install gfortran
brew install numpy
brew install scipy
```

4. Add the following to `~/.bash_profile`

```bash
export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

5. Install Python and Qt dependencies

```bash
source ~/.bash_profile
pip install numpy
pip install scipy

brew install pyqt
brew install pyqwt

pip install guiqwt
pip install guidata
```

6. Install pymzML

```bash
cd /tmp
git clone https://github.com/hroest/pymzML 
cd pymzML
sudo python setup.py install
cd ..
sudo rm -rf pymzML
```

