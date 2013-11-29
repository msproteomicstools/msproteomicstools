#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013, ETH Zurich
For a full list of authors, refer to the file AUTHORS.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""

import os

from SwathRun import SwathRun

class SwathRunCollection(object):
    """A collection of SWATH files

    Contains multiple SwathRun objects which each represent one single mass
    spectrometric injection. It can be initialized in three different ways,
    either from a set of directories (assuming each directory is one run), a
    set of files mapped to a run id (multiple files may be mapped to run id) or
    a simple flat list of chromatogram files.
    
    Attributes:
        swath_chromatograms     Dictionary mapping of the form { run_id : SwathRun}

    Methods:
        getSwathFiles:      Returns a list of SwathRun objects
        getSwathFile:       Returns a SwathRun specified by the run id
        getRunIds:          Returns all available run ids

    """

    def __init__(self):
        self.swath_chromatograms = {}

    def initialize_from_directories(self, runid_mapping):
        """Initialize from a directory
        
        This assumes that all .mzML files in the same directory are from the
        same run. There may be multiple chromatogram (chrom.mzML) files mapped
        to one run id.

        Args:
            runid_mapping(dict): a dictionary of form { run_id : directory }
        """
        self.swath_chromatograms = {}
        for runid, dname in runid_mapping.iteritems():
            import glob
            files = glob.glob(os.path.join(dname + "/*.mzML") )
            self.swath_chromatograms[ runid ] = SwathRun(files, runid)

    def initialize_from_chromatograms(self, runid_mapping):
        """Initialize from a set of mapped chromatogram files. There may be
        multiple chromatogram (chrom.mzML) files mapped to one run id.

        Args:
            runid_mapping(dict): a dictionary of form { run_id : [chromatogram_files] }
        """
        self.swath_chromatograms = {}
        for runid, chromfiles in runid_mapping.iteritems():
            self.swath_chromatograms[ runid ] = SwathRun(chromfiles, runid)

    def initialize_from_files(self, filenames):
        """Initialize from individual files, setting the runid as increasing
        integers. 

        This assumes that each .mzML file is from a separate run.

        Args:
            filenames(list): a list of filenames
        """
        self.swath_chromatograms = {}
        for i,f in enumerate(filenames):
            runid = i
            self.swath_chromatograms[ runid ] = SwathRun([f], str(runid) )

    def getSwathFiles(self):
        return self.swath_chromatograms.values()

    def getSwathFile(self, key):
        return self.swath_chromatograms[key]

    def getRunIds(self):
        return self.swath_chromatograms.keys()




