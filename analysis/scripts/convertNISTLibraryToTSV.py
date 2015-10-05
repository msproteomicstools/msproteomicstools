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

# This program takes a spectral library in NIST format and converts it to csv
#
# Input: NIST spectral library
# Output: output in csv format

import numpy as np
import sys, csv

if len(sys.argv) < 3:
    print("Please use filename as the first argument")
    exit()

print("Running NIST Parser")
fname = sys.argv[1]
outname = sys.argv[2]
fh = open(fname)

csv_headers = ['PrecursorMz', 'ProductMz', 'Tr_recalibrated', 'transition_name', 'CE',
                            'LibraryIntensity', 'transition_group_id', 'decoy', 'PeptideSequence', 'ProteinName', 
                            'Annotation', 'FullUniModPeptideName', 
                            'PrecursorCharge', 'GroupLabel', 'UniprotID', 'FragmentType', 'FragmentCharge',
                            'FragmentSeriesNumber']

writer = csv.writer(open(outname,'w'), dialect='excel-tab')
writer.writerow( csv_headers )

class Peak(object):
    def __init__(self, mz, intens, annot):
        self.mz = float(mz)
        self.intensity = float(intens)
        self.annotation = annot

rowcnt = 0
trgroup_id = 0
def handle_stack(stack):
    global rowcnt
    global trgroup_id
    Name = stack[0].split()[1]
    MW = float(stack[1].split()[1])
    comments_ = stack[2]
    npeaks = stack[3]
    comments = dict([com.split("=") for com in comments_.split() if len(com.split("=")) == 2])
    peaks = [Peak(it.split()[0], it.split()[1], it.split()[2]) for it in stack[4:] if len(it.strip()) > 0]
    print("Spectrum", Name, MW, "Nr peaks", len(peaks), npeaks)
    peaks = [p for p in peaks if p.mz > 400]
    peaks.sort(key=lambda x: x.intensity, reverse=True)
    sequence = Name.split("/")[0]
    charge = Name.split("/")[1]
    for i in range(6):
        thisrow = [comments["Parent"], peaks[i].mz, -1, rowcnt, -1, 
                  peaks[i].intensity, trgroup_id, 0, sequence, comments["Protein"], 
                  peaks[i].annotation, Name, 
                   charge, "light", comments["Protein"], -1, -1, -1
                  ]
        rowcnt += 1
        writer.writerow(thisrow)
    trgroup_id += 1

scan_cnt = 0
stack = []
started = False
MSLevel = 0
for l in fh:
    if l.find("Name") != -1:
        if started:
            handle_stack(stack)
        stack = []
        started = True
    stack.append(l)
    
handle_stack(stack)


