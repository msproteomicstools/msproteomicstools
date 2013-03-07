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
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""

import sys
import os
import csv
import getopt, glob
from configobj		import ConfigObj

from msproteomicstoolslib.data_structures.peptide 	import Peptide
from msproteomicstoolslib.format.ProteinDB			import ProteinDB  
import msproteomicstoolslib.format.speclib_db_lib	as speclib_db_lib  


def usage() :
	print ""
	print "spectrast2peakview.py"
	print ("-" * 15)
	print "This script is used as filter from spectraST files to swath input files."
	print ""
	print "Usage: "
	print "python spectrast2peakview.py [options] spectrast_file(s)"
	print "-h			--help		Display this help"
	print "-f	fasta_file	--fasta		Fasta file to relate peptides to their proteins (this is optional)."
	print "-r	iRT_file(s)	--irt	peptide iRT correspondencies"
	print "-l	mass_limits	--limits	Lower and Upper mass limits. Example: -l 400,1200"
	print "-s	ion_series	--series	List of ion series to be used. Example: -s y,b"
	print "-g	mass_modifs --gain		List of allowed fragment mass modifications. Useful for phosphorilarion. Example: -80,-98"
	print "-e			--exact		Use exact mass."
	print "-o	int		--min		Min number of reported ions per peptide/z. Default: 3"
	print "-n	int		--max		Max number of reported ions per peptide/z. Default: 20"
	print "-x	allowed_frg_charge_states	--charge	Fragment ion charge states allowed. Default: 1,2"
	print "-p   float	--precision	Maximum error allowed at the annotation of a fragment ion. Default: 0.05"
	print "-m	modifications_file	--modifications	File with the modifications delta mass"
	print "-w	swaths_file			--swaths		File containing the swath ranges. This is used to remove transitions with Q3 falling in the swath mass range."
	print "-i 	labeling_file		--isot-labeling	File containing the amino acid isotopic labeling mass shifts. If this option is used, heavy transitions will be generated."
	print "-t   time-scale			Options: minutes, seconds. Default: seconds."
	print "-d						--remove-duplicates	Remove duplicate masses from labeling"
	print ""

def writeStandardConfigFile(filename):

	config = ConfigObj()

	config.filename = filename

	# {'M': {'fv': ('-1.3094', '0'), 'sr': ('-13.0000', '0'), 'is': ('25.6100', '0'), 'pb': ('0.9500', '0'), 'ex': ('0.3476', '0'), 'hs': ('30.1900', '0'), 'sc': ('17.1900', '0')}}
	# {'S': {'pb': ('0.9774', '0'), 'xc': ('4.2610', '0'), 'dc': ('1.0000', '0'), 'fv': ('1.3848', '0')}}


	#Section SEQUEST
	config['SEQUEST'] = {}
	config['SEQUEST']['id'] = 'S'
	config['SEQUEST']['pb'] = 0.5
	config['SEQUEST']['xc'] = 2.5
	config['SEQUEST']['dc'] = 0
	config['SEQUEST']['fv'] = 0


	#Section MASCOT
	config['MASCOT'] = {}
	config['MASCOT']['id'] = 'M'
	config['MASCOT']['fv'] = 0
	config['MASCOT']['sr'] = 0
	config['MASCOT']['is'] = 0
	config['MASCOT']['pb'] = 0
	config['MASCOT']['ex'] = 0
	config['MASCOT']['hs'] = 0
	config['MASCOT']['sc'] = 0

	config.write()


def readModificationsFile(modificationsfile) :
	#Modifications is a dictionary : modifications = { modinsequence1 : ( peakview_name1 , deltaMass1 ), ... , modinsequencen : ( peakview_namen , deltaMassn )  }
	modifications = {}

	file = open(modificationsfile,"r")

	while True:
		line = file.readline()
		if len(line) == 0   : break
		if line[0]   == '#' : continue

		sline = line.split('\t')

		if len(sline) != 3 : continue

		modinsequence = sline[0]
		peakview_name = sline[1]
		deltaMass     = float(sline[2])

		modifications[modinsequence] = ( peakview_name , deltaMass )

	file.close()

	return modifications

def readiRTFile(iRT_file) :
	#Returns a dictionary of sequences and iRTs { sequence1 : iRT1 , ...}. Modifications over the sequences must be in peakview format.
	iRTs = {}

	file = open(iRT_file,"r")

	while True :
		line = file.readline()
		if len(line) == 0   : break
		if line[0]   == '#' : continue

		sequence        = ''
		iRT             = -100
		RT_experimental = 0.0

		sline = line.split('\t')

		if len(sline) == 2 :
			sequence		= sline[0]
			iRT				= float(sline[1])

		if len(sline) == 3 :
			sequence		= sline[0]
			iRT				= float(sline[1])
			RT_experimental = float(sline[2])

		iRTs[sequence] = ( iRT , RT_experimental )

	file.close()

	return iRTs

def readLabelingFile(labeling_file) :
	#Returns a dictionary of amino-acides (including also C-Term and N-Term) with the mass shifts due to an isotope labeling experiment.
	labeling = {}

	file = open (labeling_file,"r")
	while True :
		line = file.readline()
		if len(line) == 0 : break
		if line[0] == '#' : continue

		aminoacid = ''
		mass_shift = 0.0

		sline = line.split('\t')
		if len(sline) >= 2 :
			aminoacid  = sline[0]
			mass_shift = float(sline[1])

		labeling[aminoacid] = mass_shift

	file.close()

	print "Labeling file :" , labeling
	return labeling

def read_swathsfile(swathsfile) :

	swaths = []

	fsw = open( swathsfile , 'r')

	counterline = 0

	while True:
		row = fsw.readline()

		counterline += 1
		if len(row) == 0 :
			print "Swaths file contains %s swaths" % counterline
			break

		if row[0] == '#' : continue

		srow = row.split("\t")

		if len(srow) != 2 :
			print "Error when reading swaths file. Are there more than two values in the same row?"
			sys.exit(2)

		for value in srow :
			if not is_number(value) :
				print "Error when reading swaths file. Some value(s) are not numbers!"
				sys.exit(2)

		swaths.append( ( float(srow[0]) , float(srow[1]) ) )

	return swaths

def is_Q3_in_swath_range(q1 , q3 , swaths) :
	#determine which is the swath
	swath = []
	for sw in swaths :
		if q1 >= sw[0] and q1 <= sw[1] :
			swath = sw

	#If there is no swath --> return True, so that peptides out of Q1 experiment limits are removed
	if len(swath) == 0 : return True

	if ( q3 >= swath[0] and q3<=swath[1] ) : return True
	else : return False

def removeDuplicates(seq, idfun=None):
	# order preserving
	if idfun is None:
		def idfun(x): return x

	seen = {}
	result = []

	for item in seq:
		marker = idfun(item)
		if marker in seen: continue

		seen[marker] = 1
		result.append(item)

	return result


def removeSimilarDuplicates(seq, tolerance , idfun=None) :
	# order preserving
	if idfun is None:
		def idfun(x): return x

	seen = {}
	result = []

	for item in seq:
		catchup = False
		marker = idfun(item)
		if not is_number(idfun(item)) : raise "comparison values are supposed to be numbers!"

		for cp in seen.iterkeys() :
			if abs( marker - cp ) <= tolerance : catchup = True

		if catchup : continue

		seen[marker] = 1
		result.append(item)

	return result

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def filterBySearchEngineParams(searchEngineInfo, searchEngine, parameter_thresholds) :
	spectrumOK  = True
	#print parameter_thresholds
	#print searchEngineInfo

	if not 'id' in parameter_thresholds :
		return True

	if not parameter_thresholds['id'] in searchEngineInfo :
		return True

	for parameter,threshold in parameter_thresholds.iteritems() :
		if parameter in ('id') : continue
		if parameter in searchEngineInfo[parameter_thresholds['id']] :
			if float( searchEngineInfo[parameter_thresholds['id']][parameter][0] ) < float(threshold) :
				spectrumOK = False

	return spectrumOK

def translateModificationsFromSequence(sequence, modificationsLibrary) :
	'''Returns ( naked_sequence , sequence_mods_in_peakView_format , mods_in_the_Peptide_class_format ) '''
	sequence_mods_peakview = sequence

	#mods_peptide is a dictionary wich uses the position of the modification as key, and the delta mass as value:
	#example : GGGGMoxDDCDK  -> mods_peptide = { 5 : 15.99949 , 8 : 57.02147 }
	mods_peptide = {}
	mods_peptide_pv = {}
	#find (variable) modifications, expressed within brackets

	last_mod_pos = 0
	while (1):
		if '[' not in sequence[:] : break
		#print sequence
		seq_mod = ''
		#fetch the modification
		bracket0 = sequence.find('[')
		bracket1 = sequence.find(']', bracket0)
		seq_mod = sequence[bracket0+1:bracket1]
		#seq_mod has the look : 147
		if bracket0 > -1 : last_mod_pos = bracket0

		aa_modified = sequence[bracket0-1]
		aa_modified_pos = bracket0

		seq_mod = aa_modified + '[' + seq_mod + ']'
		#seq_mod has the look: M[147]

		if seq_mod in modificationsLibrary :
			pos = sequence_mods_peakview.find(seq_mod)
			end_pos = pos + len(seq_mod)
			mods_peptide[aa_modified_pos] = modificationsLibrary[seq_mod][1]
			mods_peptide_pv[aa_modified_pos] = modificationsLibrary[seq_mod][0]
		else :
			print "WARNING: Modification not known. Please, review your modifications list. This peptide modification will be ignored."
			print sequence
			#sys.exit()

		#Once processed, remove it from the sequence
		sequence_tmp = sequence[:bracket0] + sequence[bracket1+1:]
		sequence = sequence_tmp

	sequence_array = []
	for aa in sequence[:] :
		sequence_array.append(aa)

	for mod_pos,mod_txt in mods_peptide_pv.iteritems() :
		sequence_array[mod_pos-1] = sequence_array[mod_pos-1] + mod_txt

	sequence_mods_peakview = ''
	for aa in sequence_array :
		sequence_mods_peakview += aa


	sequence_no_mods =sequence

	return ( sequence_no_mods , sequence_mods_peakview , mods_peptide )

def main(argv) :

	irtfile			= ''
	fastafile		= ''
	swathsfile		= ''
	masslimits		= [0,30000]
	ionseries		= []
	useexactmass	= False
	normalization	= -1
	modificationsfile = ''
	maxtransitions	= 20
	mintransitions	= 3
	frgchargestate	= [1,2]
	precision		= 0.05
	searchEngineconfig = []
	gain_or_loss_mz_txt = []
	gain_or_loss_mz = []
	labeling = {}
	removeDuplicatesInHeavy = False
	useMinutes = False

	csv_headers = 	[	'Q1', 'Q3', 'RT_detected', 'protein_name', 'isotype',
					 'relative_intensity', 'stripped_sequence', 'modification_sequence', 'prec_z',
					 'frg_type', 'frg_z', 'frg_nr', 'iRT', 'uniprot_id', 'decoy'
					 ]

	#swaths =[	(400,426)   , (425,451)   , (450,476) , (475,501) ,
	#			(500,526)   , (525,551)   , (550,576) , (575,601) ,
	#			(600,626)   , (625,651)   , (650,676) , (675,701) ,
	#			(700,726)   , (725,751)   , (750,776) , (775,801) ,
	#			(800,826)   , (825,851)   , (850,876) , (875,901) ,
	#			(900,926)   , (925,951)   , (950,976) , (975,1001) ,
	#			(1000,1026) , (1025,1051) , (1050,1076) , (1075,1101) ,
	#			(1100,1126) , (1125,1151) , (1150,1176) , (1175,1201)
	#		]
	swaths = []

	#Get options
	try:
		opts, args = getopt.getopt(argv, "hf:l:s:en:r:m:o:w:c:z:g:i:dx:p:t:",["help","fasta","limits","series","exact","max","irt","modifications","min","swaths","config","writeconfig","gain","isot-labeling","remove-duplicates","charge","precision","timescale"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)

	argsUsed = 0
	for opt,arg in opts:
		if opt in ("-h","--help") :
			usage()
			sys.exit()
		if opt in ("-f","--fasta") :
			fastafile = arg
			argsUsed += 2
		if opt in ("-r","--irt") :
			irtfile = arg
			argsUsed += 2
		if opt in ("-m","--modifications") :
			modificationsfile = arg
			argsUsed += 2
		if opt in ("-w","--swaths") :
			print "swathsfile : " , arg
			swathsfile = arg
			argsUsed += 2
		if opt in ("-l","--limits") :
			masslimits = []
			masslimits_txt = arg.split(',')
			try :
				for val in masslimits_txt : masslimits.append( float(val) )
			except :
				print "Mass range limits are not a number! Please, try again."
				sys.exit(2)
			argsUsed += 2
		if opt in ("-s","--series") :
			ionseries = arg
			argsUsed += 2
		if opt in ("-e","--exact") :
			useexactmass = True
			argsUsed += 1
		if opt in ("-n","--max") :
			try :
				maxtransitions = int(arg)
			except :
				print "Max number of transitions is not an integer! Please, try again."
				sys.exit(2)
			argsUsed += 2
		if opt in ("-o","--min") :
			try :
				mintransitions = int(arg)
			except :
				print "Min number of transitions is not an integer! Please, try again."
			argsUsed += 2
		if opt in ("-c","--config") :
			searchEngineconfig = ConfigObj(arg)
			argsUsed+=2
		if opt in ("-z","--writeconfig") :
			writeStandardConfigFile(arg)
			sys.exit()
		if opt in ("-g","--gain") :
			gain_or_loss_mz_txt = arg.split(',')
			for obj in gain_or_loss_mz_txt :
				gain_or_loss_mz.append(float(obj))
			argsUsed += 2
		if opt in("-i","--isot-labeling") :
			labeling = readLabelingFile(arg)
			argsUsed += 2
		if opt in("-d","--remove-duplicates") :
			removeDuplicatesInHeavy = True
			argsUsed += 1
		if opt in("-x","--charge") :
			frgchargestate = []
			frgchargestate_txt = arg.split(',')
			for obj in frgchargestate_txt :
				frgchargestate.append(int(obj))
			argsUsed += 2
		if opt in("-p","--precision") :
			precision = float(arg)
			argsUsed += 2
		if opt in("-t","--timescale") :
			if arg in ["minutes","seconds"] :
				if arg in ["minutes"] : useMinutes = True
			else :
				print "Choose a right time-scale. Options are: minutes, seconds"
				sys.exit(10)



	if mintransitions > maxtransitions :
		print "This might seem a bit fool, but... You can't select a minimum number of transitions higher than the maximum!! "
		print "Min : " , mintransitions , " Max :" , maxtransitions
		sys.exit(2)



	sptxtfiles_pat = argv[argsUsed:]
	sptxtfiles = []
	for pat in sptxtfiles_pat :
		sptxtf = glob.glob(pat)
		for file in sptxtf : sptxtfiles.append(file)

	#If a modifications file is provided, read and store it into a dictionary
	modificationsLib = None
	if len(modificationsfile) > 0 :
		modificationsLib = readModificationsFile(modificationsfile)
	print "Modifications used : " , modificationsLib

	#If a fasta file is provided, read and store it into a dictionary
	print "Reading fasta file :" , fastafile
	proteins = None
	if len(fastafile) > 0 :
		proteins = ProteinDB()
		proteins.readFasta(fastafile)

	#Read iRTs file (if provided)
	iRTs = {}
	if len(irtfile) > 0 :
		iRTs = readiRTFile(irtfile)

	#Read swaths file (if provided)
	if swathsfile != '' :
		swaths = read_swathsfile(swathsfile)

	for sptxtfile in sptxtfiles :
		transitions = []
		print "Reading : " , sptxtfile
		assert sptxtfile[-6:] == '.sptxt'
		if not os.path.exists(sptxtfile):
			print "The file: %s does not exist!" % sptxtfile
			sys.exit(2)


		peakviewfilename = sptxtfile[:-6] + "_peakview.txt"
		try :
			writer = csv.writer(open(peakviewfilename,'w'), dialect='excel-tab')
		except :
			print "something went wrong while trying to write the file :" , peakviewfilename
			sys.exit(1)

		#write the headers
		writer.writerow( csv_headers )

		library_key = 99
		spectrastlib = speclib_db_lib.Library(library_key)

		num_spectrum = 0
		offset = spectrastlib.get_first_offset(sptxtfile)
		last_offset = -100

		while ( offset - last_offset > 10) :
			last_offset = offset
			offset , spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)

			#for property, value in vars(spectrum).iteritems():
			#	if property	in ['compress_spectra' ] : continue
			#	print property, ": ", value
			#sys.exit()

			sequence = spectrum.name.split('/')[0]
			z_parent = float(spectrum.name.split('/')[1])

			#print sequence, z_parent

			#find modifications in the sequence
			#mods_peptide is a dictionary wich uses the position of the modification as key, and the delta mass as value:
			#example : GGGGMoxDDCDK  -> mods_peptide = { 5 : 15.99949 , 8 : 57.02147 }
			sequence_no_mods , sequence_mods_peakview , mods_peptide  = translateModificationsFromSequence(sequence,modificationsLib)

			irt_sequence = -100
			RT_experimental = 0.0
			if spectrum.RetTime_detected != -1 :
				RT_experimental = spectrum.RetTime_detected / 60.0   #PeakView expect minutes, and spectraST reports seconds.
			if sequence_mods_peakview in iRTs :
				irt_sequence	= iRTs[sequence_mods_peakview][0]
				RT_experimental	= iRTs[sequence_mods_peakview][1]

			if useMinutes : RT_experimental = RT_experimental * 60

			spec_proteins = []
			if proteins : spec_proteins = proteins.get_proteins_containing_peptide(sequence_no_mods)


			protein_code1 = ''
			protein_desc  = ''

			for prot in spec_proteins :
				protein_code1	+= prot.code1
				protein_code1	+= ','
				protein_desc	+= prot.description
				protein_desc	+= ','


			if len(protein_code1) > 0 : protein_code1 = protein_code1[:-1]
			if len(protein_desc) > 0 :  protein_desc  = protein_desc[:-1]

			if len(protein_code1) == 0 :
				if hasattr(spectrum, 'protein_ac') : protein_code1 = spectrum.protein_ac
				else : protein_code1 = 'unknown'
			if len(protein_desc)  == 0 :
				if hasattr(spectrum, 'protein_ac') : protein_desc = spectrum.protein_ac
				else : protein_desc  = 'unknown'



			#print sequence, z_parent, protein_desc, protein_code1

			num_spectrum = num_spectrum +1

			if (num_spectrum % 1000 == 0) : print "spectra processed: %s" % num_spectrum


			pep = None
			precursorMZ = spectrum.precursorMZ
			if useexactmass : #calculate Q1 and Q3 mass/charge values
				pep = Peptide(sequence_no_mods,"",mods_peptide)
				precursorMZ = pep.getMZ(z_parent , label = '')

			#for property , value in vars(spectrum).iteritems() :
			#	print property , " : " , value

			searchenginefiltered = False
			try :
				#print spectrum.searchEngineInfo
				for searchengine in searchEngineconfig :
					if not filterBySearchEngineParams(spectrum.searchEngineInfo ,searchengine, searchEngineconfig[searchengine] ) :
						#print "filtered by search engine parameters"
						searchenginefiltered = True
						continue
			except AttributeError :
				pass



			peaks = spectrum.get_peaks()
			if searchenginefiltered : peaks = []

			filteredtransitions = []
			for peak in peaks :

				if peak.is_unknown : continue

				if peak.frg_is_isotope	: continue
				if peak.frg_z not in frgchargestate : continue
				if hasattr(peak, 'mass_error') :
					if abs(peak.mass_error) > precision : continue

				#If exact mass selected, calculate the fragment mass, otherwise keep "experimental" mass
				fragment_mz = float(peak.peak)
				if useexactmass :
					#if peak.frg_serie not in ['y','b'] :
						#for property , value in vars(peak).iteritems() :
						#	print property , " : " , value
					fragment_mz  = pep.getMZfragment(peak.frg_serie , peak.frg_nr , peak.frg_z , label = '')

				#sys.exit()


				#If ion series were specified by the user, filter those not matching user preferences.
				if len(ionseries) > 0 :
					if peak.frg_serie not in ionseries : continue

				#Filter by mass range
				if fragment_mz < masslimits[0] : continue
				if fragment_mz > masslimits[1] : continue

				#Filter fragment losses and gains, isotopes...
				if peak.is_frg_loss :
					objfound = False
					for mz in gain_or_loss_mz :
						#print mz, peak.frg_loss, abs(mz + peak.frg_loss[0] )
						if abs(mz + peak.frg_loss[0]) < 0.05 :
							objfound = True
							if useexactmass : fragment_mz += (mz / peak.frg_z )
							peak.frg_serie = peak.frg_serie + str(int(round(mz)))
					if not objfound : continue
				if peak.is_frg_gain :
					objfound = False
					for mz in gain_or_loss_mz :
						if abs(mz - peak.frg_gain[0]) < 0.05 :
							objfound = True
							if useexactmass : fragment_mz += ( mz / peak.frg_z )
							peak.frg_serie = peak.frg_serie + '+' + str(int(round(mz)))
					if not objfound : continue


				#Write the data into the data matrix (transitions)

				transition = [ precursorMZ , fragment_mz , RT_experimental , protein_desc , 'light' ,
								peak.intensity , spectrum.sequence , sequence_mods_peakview , int(z_parent) ,
								peak.frg_serie , peak.frg_z , peak.frg_nr , irt_sequence , protein_code1 , 'FALSE']

				filteredtransitions.append(transition)

			#Sort transitions by frg_serie, frg_nr, frg_z and intensity (desc) --> remove duplicates
			filteredtransitions = sorted(filteredtransitions, key= lambda x: (x[9], x[10], x[11], -x[5]))
			filteredtransitions = removeDuplicates(filteredtransitions, lambda x: (x[9], x[10], x[11]))

			#sort the transitions by intensity
			filteredtransitions = sorted ( filteredtransitions	, key=lambda transition : transition[5], reverse=True )

			#Remove transitions with very similar Q3 masses
			massTolerance = 0.02 #This should be configured by user --> TO-DO
			filteredtransitions = removeSimilarDuplicates (filteredtransitions , massTolerance , lambda x : x[1])

			#If a swaths file was provided, remove transitions falling into the swath
			filteredtransitions_tmp = []
			if len(swaths) > 0 :
				for tr in filteredtransitions :
					if not is_Q3_in_swath_range(tr[0] , tr[1] , swaths) : filteredtransitions_tmp.append(tr)
				filteredtransitions = filteredtransitions_tmp

			#if less transitions than the minimum --> continue to next spectrum
			#print filteredtransitions
			if len(filteredtransitions) < mintransitions :
				filteredtransitions = [] #I don't think this is really necessary, just in case.
				continue

			#Write in the peakview input file (until the max number of transitions per peptide/z)
			for index in range(0,min(len(filteredtransitions),maxtransitions)) :
				writer.writerow(filteredtransitions[index])

				#Isotopic labeling
				if len(labeling) > 0 :
				#heavy_transition = [ precursorMZ , fragment_mz , RT_experimental , protein_desc , 'heavy' ,
				#			peak.intensity , spectrum.sequence , sequence_mods_peakview , int(z_parent) ,
				#			peak.frg_serie , peak.frg_z , peak.frg_nr , irt_sequence , protein_code1 , 'FALSE']

					heavy_transition   	= filteredtransitions[index]
					heavy_transition[4] = 'heavy'
					precursorMZ_heavy  	= heavy_transition[0]
					fragment_mz_heavy  	= heavy_transition[1]
					sequence_heavy     	= heavy_transition[6]
					z_parent           	= heavy_transition[8]
					frg_serie          	= heavy_transition[9]
					frg_z			   	= heavy_transition[10]
					frg_number		   	= heavy_transition[11]

					#NOTE: This only works for y- and b- ions (AND their neutral losses). Other fragment series are ignored, and no heavy transition will be generated for them.
					if frg_serie[0] not in ['y','b'] : continue


					for aa in sequence_heavy :
						if aa in labeling : precursorMZ_heavy += labeling[aa] / z_parent

					frg_seq = sequence_heavy
					#b series
					if frg_serie == 'b': frg_seq = sequence_heavy[:frg_number]
					#y series
					if frg_serie == 'y': frg_seq = sequence_heavy[-frg_number:]

					for aa in frg_seq :
						if aa in labeling : fragment_mz_heavy += labeling[aa] / frg_z

					#Check for C- and N-terminal labelings
					if 'C-term' in labeling :
						precursorMZ_heavy += labeling['C-term'] / z_parent
						if frg_serie == 'y' : fragment_mz_heavy += labeling['C-term'] /frg_z
					if 'N-term' in labeling :
						precursorMZ_heavy += labeling['N-term'] / z_parent
						if frg_serie == 'b' : fragment_mz_heavy += labeling['N-term'] /frg_z


					#NOTE : if for any reason the Q3 mass has not changed (the labeling is not present in this fragment ion), it is not reported.
					#		if Q3 remains the same, then the mass shift in Q1 should be enough to fall in a different SWATH window.
					#		(In case a swaths mass ranges file is present, otherwise, we remove the not changing Q3 masses anyway).
					if removeDuplicatesInHeavy and fragment_mz_heavy == filteredtransitions[index][1] :
						if len(swaths) > 0 :
							if is_Q3_in_swath_range(precursorMZ_heavy , filteredtransitions[index][0] , swaths) : continue
						else : continue

					#Write the heavy transition into the file
					heavy_transition[0] = precursorMZ_heavy
					heavy_transition[1] = fragment_mz_heavy
					writer.writerow(heavy_transition)

			filteredtransitions = []

		print "file written : " , peakviewfilename


	print "Done."



if __name__ == '__main__':
	main(sys.argv[1:])

