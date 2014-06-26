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
$Authors: Hannes Roest, Sigurdur Smarason$
--------------------------------------------------------------------------
"""

import sys, struct, csv

"""
 *
 * Program       : Parse ABSciex Qtrap .dam files
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>, Sigurdur Smarason
 * Date          : 13.03.2012
 *
 * Note : very experimental code
          This program can parse AB Sciex Qtrap .dam files and extract
          information on the used input transition list. It will not retrieve
          any actual data, only metadata (which transitions were used).
"""

import argparse
usage  = "This program can parse AB Sciex Qtrap .dam files"
usage += "It will output a csv file that corresponds to the input transition list."

parser = argparse.ArgumentParser(description = usage )
parser.add_argument('--in', dest="inputfile", required=True, help = 'An input file (.dam format)')
parser.add_argument('--out', dest="outputfile", required=True, help = 'An input file (.csv format)')
parser.add_argument('--doAssert', action='store_true', default=False, help="Fail upon encountering an error")
args = parser.parse_args(sys.argv[1:])

inputfile = args.inputfile
outputfile = args.outputfile
do_assert = args.doAssert

class QtrapFileFormat:
    """
    Dummy class that explains the file format
    It is only based on observations and some reverse engineering. 
    There is no guarantee that this will work for newer formats of .dam files
    or, in that regard, for any .dam file.

    The Qtrap files is binary. The files have the following structure
        * Short Header, a few bytes (around 168) [part1 ]
        * At least 256 ff bytes
        * Long Header (a couple of kilobytes) [part2]
          - this header contains some information like this {{{
             Root Entry
             FileRec_StrrAcqMethodFileInfoStmxAcqMethodConfigStm&
             ()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
             AcqMethodOriginalConfigStm6
             DeviceMethod0
             DeviceMethod1
             DeviceMethod2
             DeviceMethod3
             VendorAppMethod
             VendorAppMethod
             VendorAppMethod 
             MassSpecMethodExEx&MassSpecMethodExExEx*Period0PeriodStream
             -./012345678
             Tempo LC device CH1Tempo LC device CH2Tempo nano LC Autosampler7
             D:\[...]m2.dam
             Friday March 25, 2011 03: 30: 14 
             Friday March 25, 2011 03: 30: 14 
             GenericAnalyst 1.5&File Version:  1.0
             imsbIBT-RASM11 
             D:\[...]m2.dam
             New methodTempo nano LC Autosampler
             C:\P[...]w.ini
             methodTempo LC device  CH2
             C:\P[...]w.ini
             :1New methodTempo LC device  CH1
             C:\P[...]w.ini
             :8New mass spectrometer method
             PeriodStreamEx
             Experiment0
             CycleDepParam
             ParamCollHeader
             TripleQuadMALDI
             ExperimentHeader
             ExperimentHeaderEx
             IonSourceParamsTable
             MassRangeEx
             ExperimentHeaderQTrapEX
             ExperimentHeaderQTrapEXEX4
             ExperimentHeaderDFTParamsEX8
             HpsMRM
             Parameter1
             Parameter2
             Parameter3
             Parameter4
             Parameter5
             ParameterData
             ParameterData
             ParameterData
             experiment
           }}}
        * At least 256 ff bytes
        ***********************************
        * Data / Transition list [part3]
        ***********************************
        * Padding of zero bytes
        * End of data (12x ff) with ENDOFDATASEQUENCE
        * a few hundred bytes trailing

    The data looks like this:
        * 4  bytes Q1
        * 4  bytes null
        * 4  bytes Q3
        * 4  bytes RT
        * 4  bytes null
        * 2  bytes unknown variable word (seems to be the same for each peptide)
        * xx bytes sequence in ASCII (variable, terminated by EOT)
        * 42 bytes constant bytestring (unknown )
        * 4  bytes CE (collision energy)
        * 4  bytes CE (for some reason twice, maybe for ramping?)
        * 24 bytes constant bytestring (unknown)
        == In total 96 bytes plus the length of the comment / sequence string

    """

###########################################################################
# Constants {{{
###########################################################################

BEGINOFDATASEQUENCE = "".join( [ c + '\x00' for c in 'ParameterData']) + \
        '\x00' * 38 + '\x1c\x00\x02\x01' + '\xff' * 12
        
ENDOFDATASEQUENCE = BEGINOFDATASEQUENCE   
    
Step2 = "".join( [ c + '\x00' for c in 'New experiment']) + '\xFF\xFF\x04\x00\x00\x00'
Step3 = '\x02\x00\x00\x00'+'\x04\x00\x00\x00'

EOT = '\x04' #end of transmission EOT = 04

EOTx1 = '\x04'+'\x00'
EOTx3 = '\x04'+'\x00' * 3

SEQUENCE_START_AT = 16 # counted from q1 until the sequence starts

# }}}

###########################################################################
# Start
###########################################################################

f = open(inputfile)
r = f.read()
f.close()

# Find the two header files and advance to the transition list
beginData = r.find(BEGINOFDATASEQUENCE,0)

if do_assert: assert beginData > 0

posStep3 = r.find(Step3, beginData+len(BEGINOFDATASEQUENCE)+1)

#print posStep3, beginData

if do_assert: assert posStep3>beginData


# we start in posStep3 of the file
data = r[posStep3 + 8:]
endofdata = data.find( ENDOFDATASEQUENCE )

#print endofdata

def HexStringToString(hexString):
  # convert hex string to windows 1252 string
  bytes = []
  hexStr = ''.join( hexString.split("%") )
  for i in range(0, len(hexStr), 2):
    bytes.append( chr( struct.unpack("h", hexStr[i:i+2])[0] ) )
    
  # decode as Win 1252
  string = ''.join( bytes ).decode("Windows-1252")

  return string


class Entry():

    def __init__(self): pass


class QtrapParser():

    def __init__(self, current = 0, eof = -1): 
        self.current_position = current
        self.endofdata = eof
        self.entries = []
        self.done = False

    def parse(self, data):
        entry = Entry()
        pos = self.current_position 
        
		 
		
        # 
        # Parse the data:
        entry.q1    = data[pos:pos+4]
        if data[pos:pos+4] == '\x00'*4: #if we have hit a null string here then we are done
            self.done = True; 
            return
        
        pos += 4
        null1       = data[pos:pos+4]
        pos += 4
        entry.q3    = data[pos:pos+4]
        pos += 4
        entry.rt    = data[pos:pos+4]
        pos += 4
        null2       = data[pos:pos+4]
        pos += 4+2
        # Find text (sequence) and determine break condition
        for sequence_length,c in enumerate(data[pos:]): 
            if c == EOT: break #end of transmission (end of sequence)
        
        #print sequence_length
        
        entry.sequence = data[pos:pos+sequence_length]
		
        pos += sequence_length+2 # to account for the xOO byte after it

        pos += 4  # \x44\x00\x50\x00 = DP
        entry.dp = data[pos:pos+4]
        pos += 12  # skip repeat and empyt bytes
		
        #print ":".join("{:02x}".format(ord(c)) for c in data[pos:pos+1])
        
        #make sure we are in a correct position
        if do_assert: assert data[pos:pos+2] == EOTx1
        
        pos += 2
        pos += 4  # \x45\x00\x50\x00 = EP
        entry.ep = data[pos:pos+4]
        pos += 12  # skip repeat and empyt bytes
		
        #make sure we are in a correct position
        if do_assert: assert data[pos:pos+2] == EOTx1
		
        pos += 2
        pos += 4  # \x43\x00\x45\x00 = CE
        entry.ce = data[pos:pos+4]
        pos += 12  # skip repeat and empyt bytes
		
		#make sure we are in a correct position
        if do_assert: assert data[pos:pos+2] == '\x06\x00'  #no idea why this one is different from EOT
		
        pos += 2
        pos += 6  # \x43\x00\x58\x00x50\x00 = CXP
        entry.cxp = data[pos:pos+4]
        pos += 12  # skip repeat and empty byte
        

        # sanity check whether we are at the end of the data
        if pos > self.endofdata: 
            self.done = True; 
            return
        
        self.current_position = pos  # update the current position


        # The sequence is windows 1252 encoded so we need to deal with that and convert to UTF-8
        entry.sequence = HexStringToString(entry.sequence)
        entry.sequence=entry.sequence.encode("utf-8")

        # The numbers are little-endian floats of size 4, so we use '<f' to unpack.
        entry.q1_unpacked = struct.unpack('<f', entry.q1)[0]
        entry.q3_unpacked = struct.unpack('<f', entry.q3)[0]
        entry.rt_unpacked = struct.unpack('<f', entry.rt)[0]
        entry.ce_unpacked = struct.unpack('<f', entry.ce)[0]
        entry.ep_unpacked = struct.unpack('<f', entry.ep)[0]
        entry.dp_unpacked = struct.unpack('<f', entry.dp)[0]
        entry.cxp_unpacked = struct.unpack('<f', entry.cxp)[0]
        self.entries.append(entry)


# We loop through the data part and create a new "Entry" object for each entry.
# We end the loop when we hit endofdata.
# endposition always holds the position where the new entry starts/the old one
# ends.
parser = QtrapParser(0,endofdata)
while(not parser.done):
    parser.parse(data)

# Write out a csv file
f = open(outputfile, 'wb')
csvWriter = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#write a header for the values
csvWriter.writerow( ['Q1','Q3','RT','ID','CE','EP','DP','CXP'] )
for e in parser.entries:
    csvWriter.writerow( [e.q1_unpacked, e.q3_unpacked, e.rt_unpacked, 
                        e.sequence, e.ce_unpacked,e.ep_unpacked, 
                        e.dp_unpacked, e.cxp_unpacked ] )

f.close()


# struct.unpack('<f', '\xbe\x52\xa2\x3f')
