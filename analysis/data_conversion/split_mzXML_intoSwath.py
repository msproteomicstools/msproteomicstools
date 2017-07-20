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
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""



"""
Usage: 
    python split.py filename window_size [outputdir] [noms1map]

where 
    filename is an mzXML or mzXML.gz file 
    window size is usually 32
    outputdir is . default
    noms1map is default false (ms1scans are written)

Purpose: To split files from SWATH scans into individual files.

Author: Hannes Roest loblum

###########################################################################
"""
from __future__ import print_function

import re
import sys
import os
import gzip

if len(sys.argv) < 3 or sys.argv[1] == '-h':
    print("Usage: split.py file.mzXML[.gz] windowSize [outputdir [noms1map]]]")
    sys.exit(1)

full_filename = sys.argv[1]
window_size = int(sys.argv[2])

if len(sys.argv) >= 4:
    out_dir = sys.argv[3] + "/"
    print("Set outdir to " + out_dir)
else:
        out_dir = "./"

if len(sys.argv) >= 5 and sys.argv[4] == 'noms1map':
    writeMs1 = False
    print("Skipping ms1map output")
else:
    writeMs1 = True

filename = full_filename.split(".mzXML")
if len(filename) != 2:
    print("Error, first parameter needs to be an mzXML file, was", full_filename)
    sys.exit()

filename = os.path.basename(filename[0])

# some regexes that we need since we replace parts of the XML
precursor_re = re.compile('<precursorMz([^>]*)>([^<]*)</precursorMz>')
scan_nr = re.compile('<scan num="([\d]*)"')
scan_type = re.compile('scanType="(.*)"')
ms_level = re.compile('msLevel="2"')

def rewrite_single_scan(mybuffer, swathscan):
    mybuffer = scan_nr.sub( """<scan num="%(scannr)s" """ %
            {'scannr' : swathscan}, mybuffer)
    mybuffer = scan_type.sub( """scanType="SIM" """, mybuffer)
    mybuffer = ms_level.sub( """msLevel="1" """, mybuffer)
    return mybuffer

windows = []
for w in range(window_size):
    windows.append( open(out_dir + 'split_%s_%02d.mzXML' % (filename,w) , 'w') )
if writeMs1:
    ms1map = open(out_dir + 'split_' + filename + '_' + 'ms1scan.mzXML', 'w')

# Go through all lines, find the start and end tags: <scan and </scan
# The header is found before the first <scan tag
# At each new <scan tag, a buffer is initialized to hold the current scan
# scanwindow cycles from 0 to "number of swath scans per interval"
# swathscan counts the total number of scans in each swath
scanwindow = 0
swathscan = 0
header = ''
header_done = False
mybuffer = ''
if full_filename.endswith('.gz'):
    print("Opening compressed file")
    source = gzip.open(full_filename,'rb')
else:
    source = open(full_filename)
for line in source:
    mybuffer += line
    if line.find('<scan') != -1:
        # First pass
        if not header_done:
            header_done = True
            for f in windows:
                f.write(header)
            if writeMs1:
                ms1map.write(header)
        # Start new buffer
        mybuffer = line
    if line.find('</scan') != -1:
        # We count the number of MS1 scans in the swathscan variable
        if mybuffer.find('msLevel="1"')  != -1:
            scanwindow = 0
            swathscan += 1
            if writeMs1:
                ms1map.write(mybuffer)
        else:
            mybuffer = rewrite_single_scan(mybuffer, swathscan)
            try:
                windows[scanwindow].write(mybuffer)
            except IndexError:
                print("window index %d does not fit in the array of length %d" % (scanwindow, window_size))
                raise
            scanwindow += 1
    if not header_done: header += line

for f in windows:
    f.write('  </msRun>\n</mzXML>')
    f.close()
if writeMs1:
    ms1map.write('  </msRun>\n</mzXML>')
    ms1map.close()
