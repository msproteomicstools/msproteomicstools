#!/usr/bin/env python
# -*- coding: utf-8  -*-
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
$Maintainer: Hannes Roest $
$Authors: Hannes Roest $
--------------------------------------------------------------------------
"""
from __future__ import division

# Script to generate a set of SWATH MS files for testing
# 
# This script will generate a set of files that look like SWATH MS files with
# different properties. Specifically, it will allow to se the
# IsolationWindowOffset (lower/upper) to different values to mimick what
# OpenMS/pwiz do when converting a file.
# 

import pyopenms
import numpy
from scipy.stats import norm

def getSwathExperiment(nr_scans, nr_swathes, precusorsisolation):
    exp = pyopenms.MSExperiment()
    for spec_cnt in range(1, nr_scans):
        ms1_spec = pyopenms.MSSpectrum()
        ms1_spec.setMSLevel(1)
        ms1_spec.setRT(spec_cnt*10)
        middle_scan = nr_scans // 2
        intensity = norm.pdf((spec_cnt- middle_scan) / 4.0)
        pk_list = [ [500.01, intensity*3000], [510.05, intensity*3000] ]
        peaks = numpy.array(pk_list, dtype=numpy.float32)
        ms1_spec.set_peaks(peaks)
        exp.addSpectrum(ms1_spec)
        # Swath 1: 500.01, 500.15, 500.25
        # Swath 2: 501.01, 501.15, 501.25
        # Swath 3: 502.01, 502.15, 502.25
        # Swath 4: 504.01, 504.15, 504.25
        # Swath 5: 505.01, 505.15, 505.25
        for i in range(nr_swathes):
            middle_scan = nr_scans // 2 + i # shift the middle of the gaussian by one in each scan
            intensity = norm.pdf((spec_cnt - middle_scan) / 4.0))
            intensity *= 1.2 # 20% higher intensity in each swath
            spec = pyopenms.MSSpectrum()
            spec.setMSLevel(2)
            spec.setRT(spec_cnt*10+ i+1)
            prec = pyopenms.Precursor()
            if precusorsisolation == "OpenSwath":
                prec.setIsolationWindowLowerOffset(400 + i*25)
                prec.setIsolationWindowUpperOffset(425 + i*25)
            elif precusorsisolation == "Pwiz":
                prec.setIsolationWindowLowerOffset(12.5)
                prec.setIsolationWindowUpperOffset(12.5)
            elif precusorsisolation == "Missing":
                pass
            else:
                raise Exception("precusorsisolation needs to be {Missing,Pwiz,OpenSwath}")
            prec.setMZ(400 + i *25 + 12.5);
            spec.setPrecursors( [prec])
            pk_list = [ [500.01+i, intensity*3000] , [500.15+i, intensity*3000/2.0], [500.25+i, intensity*3000/3.0]  ]
            peaks = numpy.array(pk_list, dtype=numpy.float32)
            spec.set_peaks(peaks)
            exp.addSpectrum(spec)
    return exp

exp = getSwathExperiment(20,5, "OpenSwath") 
pyopenms.MzMLFile().store("Swath_test_osw.mzML", exp)
exp = getSwathExperiment(20,5, "Pwiz") 
pyopenms.MzMLFile().store("Swath_test_pwiz.mzML", exp)
exp = getSwathExperiment(20,5, "Missing") 
pyopenms.MzMLFile().store("Swath_test_missing.mzML", exp)

