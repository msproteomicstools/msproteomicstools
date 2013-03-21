#!/usr/bin/python
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
import re, csv, sys
full_filename = sys.argv[1]
paramfilename = sys.argv[2]
outfilename = sys.argv[3]

filename = full_filename.split(".mzXML")
if len(filename) != 2:
    print "Error, first parameter needs to be an mzXML file, was", full_filename
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

def rewrite_single_scan(mybuffer, swathscan):
    # use middle point
    start = int(paramfile[ scanwindow ][0])
    end = int(paramfile[ scanwindow ][1])
    middle = (start + end)/2.0
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

