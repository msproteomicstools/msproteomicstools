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
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""
from __future__ import print_function

import sys
import os
import csv
import getopt
import operator
import random

from msproteomicstoolslib.format.ProteinDB	import	 ProteinDB  


def usage() :
	print("")
	print("proteinDBmasses.py")
	print ("-" * 15)
	print("This script retrieves the protein weights given a fasta file.")
	print("")
	print("Usage: ")
	print("python proteinDBmasses.py [options] fasta_file(s)")
	print("-h		Display this help")
	print("-a	abundancefile	Protein abundance file. If not used, the script will just report the number of peptides and molecular weight of the proteins.")
	print("-d	dynamic-range	Simulation of peptide intensity by using the peptide detectability dynamic range specified. If not used, there will not be any simulation")
	print("-e 	enzyme	Enzyme used for in-silico digestion (peptide counting). Options: trypsin, Asp-N, Arg-C, Chymotrypsin, Lys-C, Lys-N. Default: trypsin")
	print("-l 	pep-length	Minimum peptide length for the in-silico digestion. Default: 5")
	print("-m	min-peptides	Define a minimum number of peptides to identify a protein. This must be completed with the -n option.")
	print("-n	peptides-identified Define a number of identified peptides to estimate the number of proteins identified.")
	print("")


def parse_abundance(abundancefile) :
	abundances = {}
	headers = { 'protein' : "entry" , 'abundance' : "AVERAGE" }
	headers_cols = {}
	
	try: 
		rf = csv.reader(open(abundancefile,'r'), dialect='excel-tab')
	except :
		print("Something went wrong when trying to read the file : ", abundancefile)
		sys.exit(1)
	
	headerset = False
	for this_row in rf :
		if not headerset and headers['protein'] in this_row :
			headerset = True
			headers_cols = dict([ (header, idx) for idx,header in enumerate(this_row)]) 
		if not headerset : continue
		this_prot  = this_row[headers_cols[headers['protein']]]
		this_abund = this_row[headers_cols[headers['abundance']]]
		abundances[this_prot] = this_abund
	
	return abundances

def randomPepdist():
	pdf_ = [(0.003355705,1.539272326), 
		(0.003355705,2.369359294), 
		(0.003355705,3.647089191), 
		(0.013422819,5.613863463), 
		(0.030201342,8.641264671), 
		(0.030201342,13.30125957), 
		(0.036912752,20.47426076), 
		(0.046979866,31.51546298), 
		(0.046979866,48.51088001), 
		(0.053691275,74.67145511), 
		(0.040268456,114.9397044), 
		(0.036912752,176.9235061), 
		(0.104026846,272.3334568),
		 (0.097315436,419.1953535), 
		 (0.11409396,645.2558069), 
		 (0.154362416,993.2244068), 
		 (0.184563758,1528.842843)]
	
	pdf_ = [( 0.6265664160401 , 0 ),
		( 0.0012531328320802 , 1539.27232604955 ),
		( 0.0012531328320802 , 2369.35929374199 ),
		( 0.0012531328320802 , 3647.08919132535 ),
		( 0.0050125313283208 , 5613.86346284153 ),
		( 0.0112781954887218 , 8641.26467057266 ),
		( 0.0112781954887218 , 13301.2595694822 ),
		( 0.0137844611528822 , 20474.2607569056 ),
		( 0.0175438596491228 , 31515.4629794271 ),
		( 0.0175438596491228 , 48510.8800068712 ),
		( 0.0200501253132832 , 74671.4551068872 ),
		( 0.0150375939849624 , 114939.704391883 ),
		( 0.0137844611528822 , 176923.506134741 ),
		( 0.0388471177944862 , 272333.456820864 ),
		( 0.0363408521303258 , 419195.353541765 ),
		( 0.0426065162907268 , 645255.806915396 ),
		( 0.0576441102756892 , 993224.406807639 ),
		( 0.068922305764411 , 1528842.84295598 )]
	
	
	pdf =	[( 0.0102133471883268 , 0 ),
		( 0.00034935828185172 , 1539.27232604955 ),
		( 0.000537757535130571 , 2369.35929374199 ),
		( 0.000827755291951106 , 3647.08919132535 ),
		( 0.0012741408136414 , 5613.86346284153 ),
		( 0.00196124969392847 , 8641.26467057266 ),
		( 0.00301889737833723 , 13301.2595694822 ),
		( 0.00464690518965804 , 20474.2607569056 ),
		( 0.00715285256021664 , 31515.4629794271 ),
		( 0.0110101879982541 , 48510.8800068712 ),
		( 0.0169476776903155 , 74671.4551068872 ),
		( 0.0260870912595099 , 114939.704391883 ),
		( 0.0401551376428927 , 176923.506134741 ),
		( 0.0618096921224152 , 272333.456820864 ),
		( 0.0951419485656766 , 419195.353541765 ),
		( 0.146449368473575 , 645255.806915396 ),
		( 0.225425460058808 , 993224.406807639 ),
		( 0.346991172255511 , 1528842.84295598 )]
	
	cdf = [(i, sum(p for p,j in pdf if j < i)) for _,i in pdf]
	return max(i for r in [random.random()] for i,c in cdf if c <= r)

def writeMassFile(fastafile,enzyme, abundances, dynrange, minPepLength = 5, peptidesID = 0, minPepProteinID = 1) :
	base = os.path.splitext(fastafile)[0]
	massfile = base + '_mass.tsv'
	dynfile  = base + '_simulation.tsv'
	headers = ['protein_code','weight','num_of_peptides']
	headers_dyn = ['protein_code','num_of_peptides','peptide','dyn_range','protein_abundance','pep_estimated_intensity','peptide_shared_with_detected_proteins']
	peptidelib = {} # dictionary containing the simulated detectability and the protein abundances : peptides = [ detectability , { "SEQUENCE1" : { code1 : Abundance1, ...}, ... }]
	
	db = ProteinDB()
	db.readFasta(fastafile)
	try :
		writer = csv.writer(open(massfile,'w'), dialect='excel-tab')
	except :
		print("something went wrong while trying to write the file :" , massfile)
		sys.exit(1)
	if dynrange :
		try :
			writer_dyn = csv.writer(open(dynfile,'w'), dialect='excel-tab')
		except :
			print("something went wrong while trying to write the file :" , dynfile)
			sys.exit(1)
	
	if len(abundances) > 0 : headers.extend(['protein_abundance'])
	writer.writerow(headers)
	if dynrange : writer_dyn.writerow(headers_dyn)
	protein_cnt = 0
	for code, protein in db.proteinDictionary.items() :
		protein_cnt += 1
		if protein_cnt % 5000 == 0 : print("%s proteins stored" % protein_cnt)
		peptides = protein.digest(enzyme, minLength = minPepLength)
		wr_row = [protein.code2, protein.proteinWeight(), len(peptides)]
		this_abundance = -1
		if protein.code2 in abundances : 
			this_abundance = abundances[protein.code2]
			wr_row.extend([this_abundance])
		writer.writerow(wr_row)
		if dynrange and this_abundance > 0 : 
			for pep in peptides :
				if  pep not in peptidelib   : peptidelib[pep] = [ randomPepdist() , {protein : float(this_abundance) }]
				else 						: peptidelib[pep][1][protein] = float(this_abundance)
	
	
	protein_cnt = 0
	all_peptides = []
	for code, protein in db.proteinDictionary.items() :
		protein_cnt += 1
		if protein_cnt % 5000 == 0 : print("%s proteins simulated" % protein_cnt)
		if protein.code2 not in abundances : continue 
		this_abundance = -1
		this_abundance = abundances[protein.code2]
		if this_abundance < 0 : continue
		peptides = protein.digest(enzyme, minLength = minPepLength)
		for pep in peptides :
			if pep not in peptidelib : continue
			#print peptidelib[pep][0] , list(peptidelib[pep][1].itervalues())
			intensity = peptidelib[pep][0] * sum( list(peptidelib[pep][1].values()) )
			proteins_txt = ""
			if len(peptidelib[pep][1]) > 1 :
				proteins_txt = ",".join( [ str(i.code2) for i in peptidelib[pep][1].keys()] )
			wr_row = [protein.code2, len(peptides), pep, peptidelib[pep][0], this_abundance ,intensity,proteins_txt]
			all_peptides.append(wr_row)
			writer_dyn.writerow(wr_row)

	#Sort the peptides by estimated intensity
	all_peptides = sorted(all_peptides, key= lambda x: (-x[5]))
	
	protein_peptides = {}
	pep2_cnt = 0 	# Count peptides that are in proteins with at least another peptide already identified. IF it is the second peptide of the protein identified, count it twice so that we count the first one 
	pep_idx = 0
	while pep2_cnt < peptidesID and pep_idx < len(all_peptides) :
		this_peptide = all_peptides[pep_idx]
		this_sequence = this_peptide[2]
		this_protcode = this_peptide[0]
		this_intensity = this_peptide[5]
		if this_peptide[0] not in protein_peptides : protein_peptides[this_protcode] = { this_sequence : this_intensity }
		else :
			protein_peptides[this_protcode][this_sequence] = this_intensity
			if len(protein_peptides[this_protcode]) == minPepProteinID : pep2_cnt += minPepProteinID
			elif len(protein_peptides[this_protcode]) > minPepProteinID: pep2_cnt += 1
		pep_idx += 1

	#Count the number of proteins with 2 or more peptides, and proteins with only 1 peptide
	protein_2peps = 0
	protein_1pep  = 0
	for protein in protein_peptides.values() :
		#print len(protein), protein
		if len(protein) >= minPepProteinID : protein_2peps += 1
		else : 	protein_1pep += 1
	
	print("Number of proteins identified with at least %s peptides (given %s identified peptides) : %s" % (minPepProteinID,peptidesID,protein_2peps)) 
	print("Number of proteins identified with less than %s peptides (given %s identified peptides) : %s" % (minPepProteinID,peptidesID,protein_1pep)) 
	
def main(argv) :

	trypsin = {'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : ['KP', 'RP']}
	Lys_N = {'terminus' : 'N' , 'cleave' : ['K'], 'exceptions' : []}
	Asp_N = {'terminus' : 'N' , 'cleave' : ['B','D'], 'exceptions' : []}
	Arg_C = {'terminus' : 'C' , 'cleave' : ['R'], 'exceptions' : ['RP']}
	Chymotrypsin = {'terminus' : 'C' , 'cleave' : ['F','Y','W','L'], 'exceptions' : ['FP','YP','WP','LP']}
	Lys_C = {'terminus' : 'C' , 'cleave' : ['K'], 'exceptions' : ['KP']}
	
	enzymes = {'trypsin' : trypsin , 'Lys-N' : Lys_N , 'Asp-N' : Asp_N, 'Arg-C' : Arg_C , 'Chymotrypsin' : Chymotrypsin , 'Lys_C' : Lys_C}
	enzyme = trypsin
	
	abundancefile = ''
	dynrange = None
	minPepLength = 5
	#aalist = Aminoacides()
	minPepProteinID = 1
	peptidesId		= 0
	
	#Get options
	try:
		opts, _ = getopt.getopt(argv, "he:a:d:l:m:n:",["help","enzyme","abundance","dynrange","pep-length","pepprotein","pepid"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)

	argsUsed = 0
	for opt,arg in opts:
		if opt in ("-h","--help") :
			usage()
			sys.exit()
		if opt in ("-a","--abundance") :
			abundancefile = arg
			argsUsed += 2
		if opt in ("-e","--enzyme") :
			if arg not in enzymes :
				print("Error: Enzyme not recognized!") 
				print("Available enzyme options: " , [ key for key in enzymes.keys() ])
			enzyme = enzymes[arg]
			argsUsed += 2
		if opt in ("-d","dynrange") :
			dynrange = float(arg)
			argsUsed += 2
		if opt in ("-l","pep-length") :
			minPepLength = int(arg)
			argsUsed += 2
		if opt in ("-m", "pepprotein") :
			minPepProteinID = int(arg)
			argsUsed += 2
		if opt in ("-n", "pepid") :
			peptidesId = int(arg)
			argsUsed += 2

	if not os.path.exists(abundancefile) :
		print("The abundance file does not exist! " , abundancefile)
		sys.exit(2)

	abundancelib = parse_abundance(abundancefile)

	fastafiles = argv[argsUsed:]
	for fastafile in fastafiles : 
		print("processing " , fastafile)
		#File exists?
		if not os.path.exists(fastafile) :
			print("This file: %s does not exist! It will be ignored." % fastafile)
			continue
		writeMassFile(fastafile, enzyme, abundancelib, dynrange, minPepLength, peptidesID = peptidesId, minPepProteinID = minPepProteinID )
		
	pass

if __name__ == "__main__" :
	main(sys.argv[1:])