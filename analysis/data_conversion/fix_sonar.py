#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2016, ETH Zurich and Stanford University
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
This program intends to fix SONAR mzML files by adding the precursor isolation
window information into the file.

It changes the mzML by fixing the following issues:

    - replace MS-level MS:1000511 in SONAR scans with MS2 (which is wrongly written as MS1)
    - removes the MS:1000579 tag from SONAR scans (indicates MS1 scan)
    - add the precursor tag with isolation window target and offsets based on the formula provided

Usage: 
    python fix_sonar.py input.mzML paramfile output.mzML

Author: Hannes Roest

###########################################################################
"""

filename = sys.argv[1]
paramfilename = sys.argv[2]
outfilename = sys.argv[3]

if not filename.endswith("mzML"):
    print("Error, first parameter needs to be an mzML file, was", filename)
    sys.exit()

# Parameters
SONAR_CONFIG_GROUP = 2
Debug = False
WINDOW_SIZE = 22.5

## Formula to compute the center of the target window (we assume this is the
## center of the window)
def computeTarget(scanwindow):
    return 399.3983 + 2.3809 * scanwindow


# some regexes that we need since we replace parts of the XML
ms_level = re.compile('MS:1000511([^>]*)value="([^"]*)"') # MS-level
ms1_spec = re.compile('<cvParam([^>]*)MS:1000579([^>]*)>') # MS1 spec
rt_re = re.compile('MS:1000016([^>]*)value="([^"]*)"') # retention time
preset_scan_config = re.compile('preset scan configuration([^>]*)value="([\d]*)"')

outfile = open(outfilename, 'w')

# Go through all lines, find the start and end tags: <spectrum and </spectrum
# The header is found before the first <spectrum> tag
# At each new <spectrum> tag, a buffer is initialized to hold the current scan
# scanwindow cycles from 0 to "number of swath scans per interval"
# swathscan counts the total number of scans in each swath
scanwindow = 0
swathscan = 0
header = ''
header_done = False
mybuffer = ''
source = open(filename)
old_retention_time = "-1"
for line in source:
    mybuffer += line

    # Note: don't match spectrumList (!)
    if line.find('<spectrum ') != -1:
        # First pass
        if not header_done:
            header_done = True
            outfile.write(header)
        # Start new buffer
        mybuffer = line

    # End of a spectrum, write it out ...
    if line.find('</spectrum>') != -1:

        m = preset_scan_config.search(mybuffer)
        if not m:
            raise Exception("Does not look like SONAR data, no parameter 'preset scan configuration'")

        ## Determine out which kind of scan it is
        ### for Waters, group 1 is MS1
        ### for Waters, group 2 is SONAR
        ### for Waters, group 3 is lock-mass
        config_gr = int(m.group(2))
        if config_gr == SONAR_CONFIG_GROUP:

            ###################################
            # (i) Determine swath number:
            #     in SONAR, all MS2 scans have the same RT --> use this!
            m = rt_re.search(mybuffer)
            if not m:
                raise Exception("Abort, cannot find retention time")

            ###################################
            # (ii) Get retention time, compare with the one from previous scan
            retention_time = m.group(2)
            if retention_time != old_retention_time:

                if Debug:
                    print ("start new retention time scan !!!", retention_time , " != ", old_retention_time)

                # Set the new retention time, increase swathscan and re-start counter of windows
                old_retention_time = retention_time
                scanwindow = 0
                swathscan += 1


            ###################################
            # (iii) Replace the ms-level tag to use MS2 for these scans
            m = ms_level.search(mybuffer)
            if not m:
                raise Exception("Abort, cannot find MS-level")
            mybuffer = ms_level.sub("MS:1000511%svalue=\"2\"" % (m.group(1)), mybuffer)

            mybuffer = ms1_spec.sub("", mybuffer) # delete MS1 spec info

            ###################################
            # (iv) Compute precursor isolation window

            target = computeTarget(scanwindow)
            lower = WINDOW_SIZE / 2.0
            upper = WINDOW_SIZE / 2.0

            if Debug: 
                print ("we are at scanw,", scanwindow, "with swath scan", swathscan)

            # Create the precursor tag
            precursorTag = """</scanList>
  <precursorList count="1">
    <precursor>
      <isolationWindow>
        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="%(target)s" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="%(lower)s" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="%(upper)s" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
      </isolationWindow>
      <selectedIonList count="1">
        <selectedIon>
          <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="%(target)s" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
        </selectedIon>
      </selectedIonList>
    </precursor>
  </precursorList>"""  % { 'target' : target, 'lower' : lower, 'upper' : upper }
            mybuffer = mybuffer.replace("</scanList>", precursorTag)

            # Write output spectrum and increase counter
            outfile.write(mybuffer)
            scanwindow += 1
        else:
            # Its an MS1 scan or a lock-mass scan
            outfile.write(mybuffer)

    if not header_done:
        header += line

outfile.write('    </spectrumList>\n  </run>\n</mzML>\n')
outfile.close()

