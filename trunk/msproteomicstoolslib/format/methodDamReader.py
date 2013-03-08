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

import sys, struct, csv

"""
 *
 * Program       : Parse ABSciex Qtrap .dam files
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 13.03.2012
 *
 * Note : very experimental code
          This program can parse AB Sciex Qtrap .dam files and extract
          information on the used input transition list. It will not retrieve
          any actual data, only metadata (which transitions were used).
"""

from optparse import OptionParser, OptionGroup
usage = "usage: %prog inputfile outputfile [options]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "This program can parse AB Sciex Qtrap .dam files",
    "It will output a csv file that corresponds to the input transition list.")
parser.add_option_group(group)
options, args = parser.parse_args(sys.argv[1:])

inputfile = args[0]
outputfile = args[1]
do_assert = True

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

ENDOFDATASEQUENCE = "".join( [ c + '\x00' for c in 'ParameterData']) + \
        '\x00' * 38 + '\x1c\x00\x02\x01' + '\xff' * 12
AFTERSEQUENCE_CONSTANT1 = '\x04\x00D\x00P\x00\x00\x00pB\x00\x00pB' + \
        '\x00\x00\x00\x00\x04\x00E\x00P\x00\x00\x00 A\x00\x00 A' + \
        '\x00\x00\x00\x00\x04\x00C\x00E\x00'
AFTERSEQUENCE_CONSTANT2 = '\x00\x00\x00\x00\x06\x00C\x00X\x00P' + \
        '\x00\x00\x00@A\x00\x00@A\x00\x00\x00\x00'
EOT = '\x04' #end of transmission EOT = 04

STARTOFDATA = 40 # this is a magic number, where the first data structure starts
SEQUENCE_START_AT = 22 # counted from q1 until the sequence starts

# }}}

###########################################################################
# Start
###########################################################################

f = open(inputfile)
r = f.read()
f.close()

# Find the two header files and advance to the transition list
part1 = r.find( '\xff' * 256, 0)
part2 = r.find( '\xff' * 256, part1 + 256)
part3 = r.find( '\xff' * 256, part2 + 256)
if do_assert: assert r.find( '\xff' * 256, part3 + 256) == -1

# Find the end of the \xff sequence (and the start of the transition list)
for start,p in enumerate(r[part3:]): 
    if p != '\xff': break

# we start in part 3 of the file after the block of ff's 
data = r[part3 + start:]
endofdata = data.find( ENDOFDATASEQUENCE )

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
        sequence_start = self.current_position + SEQUENCE_START_AT 

        # 
        # Parse the data:
        #   Q1, null, Q3, RT, null, null
        #   text
        #   CE (collision energy)
        # 
        entry.q1    = data[pos:pos+4]
        null1       = data[pos+4:pos+8]
        entry.q3    = data[pos+8:pos+12]
        entry.rt    = data[pos+12:pos+16]
        null2       = data[pos+16:pos+20]
        entry.unknown_variable_word = data[pos+20:pos+22]

        # Find text (sequence) and determine break condition
        for sequence_length,c in enumerate(data[sequence_start:]): 
            if c == EOT: break #end of transmission (end of sequence)
        entry.sequence = data[sequence_start:sequence_start+sequence_length]

        # Calculate the position of the collision energy, store the unknown 42
        # bytes in afters1 and the unknown 24 bytes in afters3. Store the
        # collision energy (CE) and advance the current position.
        ceposition = sequence_start + sequence_length + 42 
        entry.afters1 = data[sequence_start+sequence_length:ceposition]
        entry.ce = data[ceposition:ceposition+4]
        entry.ce_repeat = data[ceposition+4:ceposition+4+4]
        entry.afters3 = data[ceposition+4+4:ceposition+32]
        self.current_position = ceposition + 32  

        # check whether we are done parsing all entries
        if ceposition > self.endofdata: 
            self.done = True; 
            return

        # Assert that some of our assumptions about the file are true
        if do_assert: assert null1 == '\x00' * 4
        if do_assert: assert null2 == '\x00' * 4
        if do_assert: assert entry.afters1 == AFTERSEQUENCE_CONSTANT1
        if do_assert: assert entry.ce == entry.ce_repeat
        if do_assert: assert entry.afters3 ==  AFTERSEQUENCE_CONSTANT2

        # The sequence is spaced with zero bytes, so we just remove them.
        # The numbers are little-endian floats of size 4, so we use '<f' to unpack.
        entry.sequence = "".join( [s for s in entry.sequence if s != '\x00'] )
        entry.q1_unpacked = struct.unpack('<f', entry.q1)[0]
        entry.q3_unpacked = struct.unpack('<f', entry.q3)[0]
        entry.rt_unpacked = struct.unpack('<f', entry.rt)[0]
        entry.ce_unpacked = struct.unpack('<f', entry.ce)[0]
        self.entries.append(entry)

# We loop through the data part and create a new "Entry" object for each entry.
# We end the loop when we hit endofdata.
# endposition always holds the position where the new entry starts/the old one
# ends.
parser = QtrapParser(STARTOFDATA, endofdata)
while(not parser.done):
    parser.parse(data)

# Write out a csv file
f = open(outputfile, 'wb')
csvWriter = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
for e in parser.entries:
    csvWriter.writerow( [e.q1_unpacked, e.q3_unpacked, e.rt_unpacked, 
                        e.sequence, e.ce_unpacked] )

f.close()


# struct.unpack('<f', '\xbe\x52\xa2\x3f')
