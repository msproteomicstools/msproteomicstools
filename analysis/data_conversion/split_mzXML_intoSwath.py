#!/usr/bin/env python
# -*- coding: utf-8  -*-

# // -*- mode: C++; tab-width: 2; -*-
# // vi: set ts=2 fdm=marker:

"""
 *
 * Program       : SWATH - split.py
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 23.05.2012 
 *
 *
 * Copyright (C) 2012 Hannes Roest
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
 *
"""

"""
Usage: 
    python split.py filename window_size [outputdir] [noms1map]

where 
    filename is an mzXML file 
    window size is usually 32
    outputdir is . default
    noms1map is default false (ms1scans are written)

Purpose: To split files from SWATH scans into individual files.

Author: Hannes Roest loblum

###########################################################################
"""

import re, csv, sys, os
full_filename = sys.argv[1]
window_size = int(sys.argv[2])

if len(sys.argv) >= 4:
    out_dir = sys.argv[3] + "/"
    print "Set outdir to " + out_dir
else:
        out_dir = "./"

if len(sys.argv) >= 5 and sys.argv[4] == 'noms1map':
    writeMs1 = False
    print "Skipping ms1map output"
else:
    writeMs1 = True

filename = full_filename.split(".mzXML")
if len(filename) != 2:
    print "Error, first parameter needs to be an mzXML file, was", full_filename
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
    windows.append( open(out_dir + 'split_' + filename + '_' + str(w) + '.mzXML', 'w') )
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
            windows[scanwindow].write(mybuffer)
            scanwindow += 1
    if not header_done: header += line

for f in windows:
    f.write('  </msRun>\n</mzXML>')
    f.close()
if writeMs1:
    ms1map.write('  </msRun>\n</mzXML>')
    ms1map.close()

