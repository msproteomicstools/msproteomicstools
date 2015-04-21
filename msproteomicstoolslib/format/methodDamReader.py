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

class QtrapFileFormat:
    """
    Dummy class that explains the file format
    It is only based on observations and some reverse engineering. 
    There is no guarantee that this will work for newer formats of .dam files
    or, in that regard, for any .dam file.

    The Qtrap files is binary. Entries are in a multiple of 128 bytes and
    all start with a name, i.e. a sequence of ascii characters separated by a null
    byte except the header at the beginning which has no name but is 1024 bytes
    long. Entries vary with the equipment profile (e.g. type of LC) as well as
    the experimental profile. If entries surpass 128 bytes or the whole file is
    not a multiple of 128 bytes it is padded with a series of xFF bytes in various
    places. Some of these paddings can occur within entries. The begin and end of
    these paddings are not always clear although many times the end is signaled with 
    x04x00. xFF seems also to be used as a delimiter and/or as mini padding in
    entries.
    
    The following are named entries that have been found, some entries are 
    found multiple times and there seems to be some sequence to the entries
    but right now this sequence is not fully known.
    
    	Root Entry {two per file}
		FileRec_Str
		AcqMethodFileInfoStm
		AcqMethodConfigStm
		AcqMethodOriginalConfigStm
		DeviceMethod0
		DeviceMethod1  {possibly more if there are more devices in the instrument profile?}
		LCPumpMethod { specific to one vendor? }
		ShimadzuMethod { specific to one vendor }
		ShimadzuVendorMethod { specific to one vendor }
		ShimadzuPretreatment { specific to one vendor }
		MSConfigInfo
		MSConfigInfoExp
		MSConfigInfoDMS
		MSConfigInfoDualMass
		MassSpecMethod
		MassSpecMethodEx
		MassSpecMethodExEx
		MassSpecMethodDMS
		MassSpecMethodExExEx
		Period0
		Period1 {possibly more if there are more periods in the method?}
		{
		Controller
		Pumps 
		TimeProgram
		Detectors
		Autosampler
		PDA
		OVEN
		SUBCONTROLLER
		Pretreatment
		Text
		} These might be vendor and/or equipment specific
		PeriodStream
		PeriodStreamEx
		Experiment0
		Experiment1 {possibly more if there are more experiments in the method?}
		CycleDepParam
		ParamCollHeader
		{
		TripleQuadMALDI
		ExperimentHeader
		ExperimentHeaderEx
		IonSourceParamsTable
		MassRangeEx
		MassRangeExEx
		ExperimentHeaderQTrapEX
		ExperimentHeaderQTrapEXEX
		ExperimentHeaderDFTParamsEX
		ExperimentHeaderFJ
		sMRM
		sMRMEX
		} possibly more sets of these entries depending on the method?
		{
		MassRangeExEx
		MassRangeEx
		ParamCollHeader
		Parameter0
		Parameter1
		Parameter2
		Parameter3
		Parameter4
		Parameter5
		}  possibly more sets of these entries depending on the method?
		{
		ParameterData
		ParameterData
		ParameterData
		ParameterData
		ParameterData
		ParameterData
		}  possibly more sets of these entries depending on the method?
		
		MRM data can be found in one or more of the Parameter entries
		and is indicated with the x02x00x00x00x04x00x00x00 sequence. 
		Immediately after that the MRM data starts. In general the data
		looks like this
        * 4  bytes Q1
        * 4  bytes null
        * 4  bytes Q3
        * 4  bytes RT (in min) or dwell time (in ms) depending on the method
        * 4  bytes null
        * 2  bytes unknown
        * xx bytes ID probably in Windows 1252 encoded text (variable, terminated by x04x00) with x00 between
        * the characters
        * DP,EP,CE,CXP will each be spelled out in ascii (with the usual x00 between) followed by the
        begin and end values (four bytes for each parameter). Possibly these can be use to ramp the settings although in the method
        setup on the Qtrap it is not immediately obvious how to do so. Note that padding of xFF can occur between these parameters.
        The next Q1 value if present will be found 6 (length of CXP) + 12 bytes after the position of the CXP string.
        If there isn't one then the value of the Q1 bytes in that position will be x00x00x00x00
        
        Other source parameters such as curtain gas (CUR), source temperature (TEM), GS1, GS2, CAD and IS
        can also be found spelled out in the paramter entries. Sometimes in the same entry that contains the MRM 
        data and sometimes not.  They also seem to have start and stop values (4 bytes each) and are sometimes 
        repated more than once in the same paramter entry and/or between parameter entries.

    """

###########################################################################
# Constants {{{
###########################################################################


if sys.version_info >= (3,0,0):
    ParameterDatafield = b"".join( [ bytes(c, "ascii") + b'\x00' for c in 'ParameterData']) + \
            b'\x00' * 38 + b'\x1c\x00\x02\x01' + b'\xff' * 12

    Q1start = b'\x02\x00\x00\x00'+b'\x04\x00\x00\x00'

    EOT = b'\x04' #end of transmission EOT = 04
else:
    ParameterDatafield = "".join( [ c + '\x00' for c in 'ParameterData']) + \
            '\x00' * 38 + '\x1c\x00\x02\x01' + '\xff' * 12
                
    Q1start = '\x02\x00\x00\x00'+'\x04\x00\x00\x00'
     
    EOT = '\x04' #end of transmission EOT = 04

# }}}

###########################################################################
# Start
###########################################################################


def HexStringToString(hexString):
	# convert hex string to windows 1252 string
	bytes = []
	hexStr = ''.join( hexString.split("%") )
	
	#these can not be converted to characters
	badCodes = [127, 129, 141, 143, 144, 157]
	
	for i in range(0, len(hexStr), 2):
		if len(hexStr[i:i+2])==2:
			dum=struct.unpack("<H", hexStr[i:i+2])[0]
		else:
			dum=-1

		test = True
		if dum<32:
			test = False
			
		if dum>255:
			test = False
		
		if dum in badCodes:
		 	test = False
		 	
		
		if test == True:
			bytes.append( chr( dum ) )
		else:
			#special check for the xFF value in order to better understand the
			#layout of the data
			if dum == 65535:
				bytes.append( chr ( 70) )
			else:
				bytes.append( chr (46) )
		
    
	# decode as Win 1252
	string = ''.join( bytes ).decode("Windows-1252")

	return string



class Entry():

    def __init__(self): pass


class MRMParser():

	def __init__(self, data): 
		self.current_position = 0
		self.data = data
		self.endofdata = len(data)
		self.entries = []
		self.done = False
		self.ionArr = ['DP','EP','CE','CXP']	
		self.ionParams = {'DP':'\x00'*4,'EP':'\x00'*4,'CE':'\x00'*4,'CXP':'\x00'*4}
	
	def parse(self):
		self.current_position = self.data.find(Q1start)+len(Q1start)
		if self.current_position>-1:
			while(not self.done):
				self.parseOne()
				
	
	def parseOne(self):
		entry = Entry()
		pos = self.current_position 
	
		# 
		# Parse the data:
		entry.q1    = self.data[pos:pos+4]
		if entry.q1 == '\x00'*4: #if we have hit a null string here then we are done
			self.done = True; 
			return
        
		pos += 4 #move 4 positions ahead
		pos += 4 # 4x null byte
		entry.q3    = self.data[pos:pos+4]
		pos += 4 #move 4 positions ahead
		entry.rt    = self.data[pos:pos+4]
		pos += 4 #move 4 positions ahead
		pos += 4 # 4x null byte
		pos += 2 # 2x unknown byte
		# Find text (sequence) and determine break condition
		for sequence_length,c in enumerate(self.data[pos:]): 
			if c == EOT: break #end of transmission (end of sequence)
        
		#print sequence_length
        
		entry.sequence = self.data[pos:pos+sequence_length]
		
		pos += sequence_length+2 # to account for the xOO byte after it
        
		#look for the DP, EP, CE and CXP, data structure splits on CXP+12 bytes before next Q1
		ionParam=self.ionParams
		for e in self.ionArr:
			if sys.version_info >= (3,0,0):
				myStr = b"".join( [ bytes(c, "ascii") + b'\x00' for c in e])
			else:
				myStr = "".join( [ c + '\x00' for c in e])
			foundOne = self.data[pos:].find(myStr)
			if foundOne>-1:
				pos += foundOne + len(myStr)
				ionParam[e] = self.data[pos:pos+4]
		if foundOne>-1:  #CXP was found...
			pos += 12

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
		entry.ce_unpacked = struct.unpack('<f', ionParam[b'CE'])[0]
		entry.ep_unpacked = struct.unpack('<f', ionParam['EP'])[0]
		entry.dp_unpacked = struct.unpack('<f', ionParam['DP'])[0]
		entry.cxp_unpacked = struct.unpack('<f', ionParam['CXP'])[0]
		self.entries.append(entry)

def main(args):

    inputfile = args.inputfile
    outputfile = args.outputfile
    do_assert = args.doAssert

    f = open(inputfile, 'rb')
    r = f.read()
    f.close()

    # Find all instances of ParameterDatafield and collect the resulting data
    beginData = r.find(ParameterDatafield, 0)
    endData = beginData
    if do_assert: assert beginData > 0

    #loop to put all ParameterData entries into an array
    paramData = []
    end=False
    while(endData>-1):
        endData = r.find(ParameterDatafield,beginData+len(ParameterDatafield)+1)
        if endData>beginData:
            paramData.append(r[beginData:endData])
        elif endData==-1:
            paramData.append(r[beginData:])
        beginData=endData
        

    #look for the begin tag of the Q1/Q3/time packet in each Parameterdata entry and parse it for MRM transitions if found
    MRMs = []
    for e in paramData:
        beginData=e.find(Q1start)
        if beginData>-1:
            parser = MRMParser(e)
            parser.parse()
            MRMs.append(parser.entries)

        
    # Write out a csv file
    if sys.version_info >= (3,0,0):
        f = open(outputfile, 'w', newline='')
    else:
        f = open(outputfile, 'wb')

    csvWriter = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    #write a header for the values
    csvWriter.writerow( ['Q1','Q3','Time','ID','CE','EP','DP','CXP'] )
    for oneMRMpacket in MRMs:
        for e in oneMRMpacket:
            csvWriter.writerow( [e.q1_unpacked, e.q3_unpacked, e.rt_unpacked, 
                    e.sequence, e.ce_unpacked,e.ep_unpacked, 
                    e.dp_unpacked, e.cxp_unpacked ] )

    f.close()

def handle_args():
    import argparse
    usage  = "This program can parse AB Sciex Qtrap .dam files"
    usage += "It will output a csv file that corresponds to the input transition list."

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="inputfile", required=True, help = 'An input file (.dam format)')
    parser.add_argument('--out', dest="outputfile", required=True, help = 'An input file (.csv format)')
    parser.add_argument('--doAssert', action='store_true', default=False, help="Fail upon encountering an error")
    args = parser.parse_args(sys.argv[1:])
    return args

if __name__=="__main__":
    options = handle_args()
    main(options)
