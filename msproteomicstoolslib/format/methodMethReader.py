#!/usr/bin/python
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

import sys, StringIO, csv


"""
 *
 * Program       : Parse Thermo TSQ .meth files
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 13.03.2012
 *
"""

from optparse import OptionParser, OptionGroup
usage = "usage: %prog inputfile outputfile [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "This program can parse Thermo TSQ .meth files",
    "It will output a csv file that corresponds to the input transition list.")
parser.add_option_group(group)
options, args = parser.parse_args(sys.argv[1:])

inputfile = args[0]
outputfile = args[1]

class TSQFileFormat:
    """
    The TSQ files is mixed binary and text. The files have the following structure
        * Header [part1 ]
        * At least 256 ff bytes
        * Header [part2 ]
        * At least 256 ff bytes
        * Data part, looks like this:
    TSQ Quantum Instrument Method

    Creator: Quantum Last modified: 3/24/2011 08:06:10 AM by Quantum


    TSQ MS Method Settings:

    Method Type: EZ Method
    MS Run Time (min): 70.00
    Experiment Type: SRM
    Chrom Filter Peak Width (s): Not used
    Collision Gas Pressure (mTorr): 1.5
    Use Tuned Tube Lens Value: Yes
    Q1 Peak Width (FWHM): 0.70
    Display Time Range for SRM table: Yes
    Cycle Time (s): 2.000
    Skimmer Offset (V): Not used

    SRM Table:
     Parent Product CE Start Stop Pol Name
     482.273 589.334 25 0.00 50.00 + FGLEGLESVVPGIK
     [ ... data ...] 

    Tune Method - Adjustable parameters:
     Capillary Temperature: 280.0
     Vaporizer Temperature: 0.0
     Sheath Gas Pressure: 0.0
     Ion Sweep Gas Pressure: 0.0
     Aux Valve Flow: 0.0
     Spray Voltage: Positvie palority - 1300.0 , Negative polarity - 3000.0
     Discharge Current: Positvie palority - 4.0 , Negative polarity - 4.0

    Divert Valve: not used during run
    """

f = open(inputfile, 'rb')
r = f.read()
f.close()

# there are three parts of the data format, separated by a block of 
# 256 (or more) ffs
part1 = r.find( '\xff' * 256, 0)
part2 = r.find( '\xff' * 256, part1 + 256)
part3 = r.find( '\xff' * 256, part2 + 256)
assert r.find( '\xff' * 256, part3 + 256) == -1

# find the start of the data
for start,p in enumerate(r[part3:]): 
    if p != '\xff': break

# get all data and remove zero bytes and extra spaces
data = r[part3+start:] 
cleandata = "".join( [p for p in data if p != '\x00'] )
oldlen = len(cleandata) + 1
while(oldlen != len(cleandata)):
    oldlen = len(cleandata)
    cleandata = cleandata.replace( '  ', ' ')

# read data as csv file and store again to write as "real" csv file
source = StringIO.StringIO( cleandata )
reader = csv.reader(source, delimiter=' ')
rows = list(reader)
output = []
assert rows[19] == ['', 'Parent', 'Product', 'CE', 'Start', 'Stop', 'Pol', 'Name']
for row in rows[20:]: 
    if len(row) == 0: break
    output.append([ row[1], row[2], ( float(row[4]) + float(row[5]) ) /2.0,
        row[7], row[3] ])

f = open(outputfile, 'wb')
csvWriter = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
for o in output: csvWriter.writerow(o)
f.close()

