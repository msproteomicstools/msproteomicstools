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
from __future__ import division
from __future__ import print_function

import re
import csv
import sys

"""
This program intends to fix SWATH mzXML files by adding the precursor isolation
window information into the file.

Usage: 
    python split.py filename paramfile

filename is an mzXML file. The paramfile is tab-delimited file that contains
the start and end m/z in the 1st and 2nd column and has as many rows as there
are windows.

Author: Hannes Roest

###########################################################################
###########################################################################
Paramfile Format
{{{ 
400	427	2
423	452	3
...

}}}
###########################################################################
###########################################################################
"""

# The paramfile contains the swathes and is needed here to get the correct
# number of swath-windows.
full_filename = sys.argv[1]
paramfilename = sys.argv[2]
outfilename = sys.argv[3]

filename = full_filename.split(".mzXML")
if len(filename) != 2:
    print("Error, first parameter needs to be an mzXML file, was", full_filename)
    sys.exit()

filename = filename[0]

# parse parameter file as csv tab delimited
f = open(paramfilename)
r = csv.reader(f, delimiter='\t')
paramfile = [l for l in r]
f.close()

# some regexes that we need since we replace parts of the XML
precursor_re = re.compile('<precursorMz([^>]*)>([^<]*)</precursorMz>')
scan_nr = re.compile('<scan num="([\d]*)"')
scan_type = re.compile('scanType="(.*)"')
ms_level = re.compile('msLevel="2"')

def estimateCorrectSwathScan(center):
    # Try to estimate which was the correct Swath window from the center of the
    # window alone.

    result = []
    for entry in paramfile:
        if center >= float(entry[0]) and center <= float(entry[1]):
            result.append([float(entry[0]), float(entry[1])])
    
    if len(result) > 2:
        raise Exception("More than one entry in the param file match the center", center)
    if len(result) == 0:
        raise Exception("No entry in the param file matches the center", center)
    
    return result[0]

def rewrite_single_scan(mybuffer, swathscan):
    # use middle point
    try:
        start = float(paramfile[scanwindow][0])
        end = float(paramfile[scanwindow][1])
        middle = (start + end) / 2.0
        width = end - start
    except IndexError:
        # Catch error condition
        start = -1
        end = -1

    # TODO check if window wideness is already present ... 
    precursor_match = precursor_re.search(mybuffer)
    if precursor_match:
        # Ensure that the center of the window is inside the SWATH 
        old_center = float(precursor_match.group(2))
        if not (old_center > start and old_center < end):
            print("WARNING: previous precursorMz", old_center, \
              "did not fall inside the expected SWATH window %s to %s." % (start, end))
            start, end = estimateCorrectSwathScan(old_center)
            print("  Will replace it with a window from %s to %s" % (start, end))
            middle = (start + end) / 2.0
            width = end - start

    mybuffer = precursor_re.sub( """<precursorMz windowWideness="%(width)s"\\1>%(middle)s</precursorMz>""" %
            {'width' : width, 'middle' : middle }, mybuffer)
    mybuffer = scan_type.sub( """scanType="SIM" """, mybuffer)
    ## Removed the renumbering of the scans and the replacement of MS2 with MS1 level
    ## mybuffer = scan_nr.sub( """<scan num="%(scannr)s" """ %
    ##         {'scannr' : swathscan}, mybuffer)
    ## mybuffer = ms_level.sub( """msLevel="1" """, mybuffer)
    return mybuffer

outfile = open(outfilename, 'w')

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
source = open(filename + '.mzXML')
for line in source:
    mybuffer += line
    if line.find('<scan') != -1:
        # First pass
        if not header_done:
            header_done = True
            outfile.write(header)
        # Start new buffer
        mybuffer = line
    if line.find('</scan') != -1:
        # We count the number of MS1 scans in the swathscan variable
        if mybuffer.find('msLevel="1"')  != -1:
            scanwindow = 0
            swathscan += 1
            outfile.write(mybuffer)
        else:
            mybuffer = rewrite_single_scan(mybuffer, swathscan)
            outfile.write(mybuffer)
            scanwindow += 1
    if not header_done: header += line

outfile.write('  </msRun>\n</mzXML>')
outfile.close()
