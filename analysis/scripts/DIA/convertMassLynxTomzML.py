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
from __future__ import print_function

import pyopenms
import numpy as np
import sys

# This program takes mass spectrometric data in ASCII format as produced by
# Waters DataBridge program (V4.1) with the MassLynx to ASCII option and
# convets it to mzML
#
# Input: ASCII data
# Output: output in mzML format

if len(sys.argv) < 3:
    print("Please use filename as the first argument and the outputname as the second argument")
    exit()

fname = sys.argv[1]
outname = sys.argv[2]
fh = open(fname)
outexp = pyopenms.MSExperiment()

def handle_stack(stack, mslevel):
    # print "Size", len(stack), mslevel
    scan_nr_ = stack[0].split()
    assert(len(scan_nr_) == 2)
    scan_nr = int(scan_nr_[1])
    rt_ = stack[1].split()
    assert(len(rt_) == 3)
    rt = float(rt_[2])*60
    # Convert to mz/int pairs
    pairs = [it.split() for it in stack[3:] if len(it.strip()) > 0 and float(it.split()[1]) > 0.0]
    #mz = [float(it.split()[0]) for it in stack[3:] if len(it.strip()) > 0]
    #intensity = [float(it.split()[1]) for it in stack[3:] if len(it.strip()) > 0]
    try:
        mz = [float(it[0]) for it in pairs]
        intensity = [float(it[1]) for it in pairs]
    except ValueError:
        print("Could not convert", len(stack), "with pairs" , len(pairs))
        return
    assert len(mz) == len(intensity)
    #
    print("Handle scan at rt", rt)
    peaks = np.ndarray(shape=(len(mz), 2), dtype=np.float32)
    peaks[:,0] = mz
    peaks[:,1] = intensity
    s = pyopenms.MSSpectrum()
    s.set_peaks(peaks)
    s.setRT(rt)
    s.setMSLevel(1)
    outexp.addSpectrum(s)

scan_cnt = 0
stack = []
started = False
MSLevel = 0
for l in fh:
    if l.find("Scan") != -1:
        if started:
            #break
            handle_stack(stack, MSLevel+1)
            MSLevel = (MSLevel+1)%2
        stack = []
        started = True
    stack.append(l)
    
pyopenms.MzMLFile().store(outname, outexp)


