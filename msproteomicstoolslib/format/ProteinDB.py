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

from __future__ import print_function
import os
import csv
import sys
import getopt
import re
try:
    # Python 2
    from itertools import izip
except ImportError:
    # Python 3
    izip = zip

from ..data_structures import aminoacides

class Protein :
    """Represents a protein which may be stored in a ProteinDB instance"""
    
    def __init__(self , proteindetails = None, calWeight = False) :
        self.code1 = ''
        self.code2 = ''
        self.modres = ''
        self.ncbi_tax_id = ''
        self.description = ''
        self.sequence = ''
        
        if proteindetails :
            if 'code1'            in proteindetails : self.code1            = proteindetails['code1']
            if 'code2'            in proteindetails : self.code2            = proteindetails['code2']
            if 'modres'            in proteindetails : self.modres            = proteindetails['modres']
            if 'ncbi_tax_id'    in proteindetails : self.ncbi_tax_id    = proteindetails['ncbi_tax_id']
            if 'description'    in proteindetails : self.description    = proteindetails['description']
            if 'sequence'        in proteindetails : self.sequence        = proteindetails['sequence']
        
        if calWeight : self.weight = self.proteinWeight()

    def digest(self, cleavageRules = {'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : ['KP', 'RP']}, minLength = 0, missedCleavages = 0):
        peptides = []
        run_seq = 1
        if cleavageRules['terminus'] == 'N' : run_seq = -1
        
        grouped = lambda iterable,n : izip(*[iter(iterable)]*n)
        
        #Build the regular expression
        cleaveregex = '' 
        for cleave in cleavageRules['cleave'] :
            this_regex = '(' + cleave + ')'
            for excep in cleavageRules['exceptions'] :
                if cleave in excep :
                    if excep.find(cleave) == 0 : this_regex = this_regex + '(?!' + excep[1:] + ')'  #It is a C-terminus exception (i.e. a proline AFTER the cleavage site)
                    else : this_regex = '(?!' + excep[:-1] + ')' + this_regex #It is an N-terminus exception
                    excep_cleave = cleavageRules['exceptions']
            if len(cleaveregex) > 0 : cleaveregex = cleaveregex + '|' 
            cleaveregex = cleaveregex + this_regex

        mm = re.split(cleaveregex, self.sequence)
        while None in mm : mm.remove(None)
        
        for x,y in grouped(mm[::run_seq],2) :
            this_peptide = ''
            if run_seq == 1  : this_peptide = x+y
            if run_seq == -1 : this_peptide = y+x  
            peptides.append( this_peptide )
        
        pos_remaining = len(mm) - 1 
        if run_seq == -1 : pos_remaining = 0
        if mm[pos_remaining] not in cleavageRules['cleave'] : peptides.append( mm[pos_remaining] )

        mCleav = []
        for i in range(1,missedCleavages+1) :
            for j,pep in enumerate(peptides) :
                if j < len(peptides) - i : 
                    curr_missedCleav = ''.join(peptides[j:j+i+1])
                    mCleav.append(curr_missedCleav) 

        #Remove peptides smaller than minimum
        final_peptides = []
        for p in peptides : 
            if len(p) >= minLength : final_peptides.append(p) 
        for p in mCleav :
            if len(p) >= minLength : final_peptides.append(p)

        return final_peptides
        
        
    def proteinWeight(self) :
        weight = 0.0
        aaLib = aminoacides.Aminoacides()
        for aa1 in self.sequence[:] :
            for aa2 in aaLib.list :
                if aa2.code == aa1 : 
                    weight += aa2.deltaMass
                    break
        return weight
    
    def pseudoreverse(self,decoytag = "DECOY_") :

        cleaveage_aa = ['K','R']
    
        index = 0
        rev_peptide = ''
        rev_protein = ''
        current_peptide = ''
        while True :
            if index > len(self.sequence) - 1 : break  #End of the protein
            if self.sequence[index] in cleaveage_aa : 
                #reverse the current peptide and add it to the reversed protein
                for aa in current_peptide[::-1] :
                    rev_peptide += aa
                rev_peptide += self.sequence[index]
                rev_protein += rev_peptide
                rev_peptide = ''    
                current_peptide = ''
            else :
                current_peptide += self.sequence[index]
            index+=1
        
        #Reverse the last tail
        if len(rev_protein) < len(self.sequence) :
            rev_protein += self.sequence[:len(rev_protein)-1:-1]
        
        #If there is not any cleavage site at the protein --> reverse the total protein
        if len(rev_protein) == 0 : rev_protein = self.sequence[::-1]
        
        if not len(rev_protein) == len(self.sequence) :
            print("Warning : size of reversed protein does not match with the original. " , len(rev_protein) , len(self.sequence) )
                
        self.sequence = rev_protein
        self.code1    = decoytag + self.code1

class ProteinDB :
    """Represents a set of proteins."""
    
    def __init__(self, dict=None) :
        if dict : self.proteinDictionary = dict
        else : self.proteinDictionary = {}
        self.chunksize = 80
    
    def findSequenceInProteins(self, sequence) :
        '''Given a peptide sequence, this def returns a list of proteins (code1) in which the peptide is present. Don't use modified sequences, remove their modifications before calling this def.'''
        proteinList = []
        
        for protein_code1 , protein in self.proteinDictionary.iteritems() :
            if sequence in protein.sequence : proteinList.append(protein_code1)
    
        return proteinList

    def pseudoreverseDB(self,decoytag = "DECOY_") :
    
        for protein_code1, protein in self.proteinDictionary.iteritems() :
            protein.pseudoreverse(decoytag)

    def writeFastaFile(self, fastaFileName, chunksize = -1, format=None) :
        if chunksize < 0 : chunksize = self.chunksize - 1
        file = open(fastaFileName,"w")
        
        for code1, protein in self.proteinDictionary.iteritems() :
            if not format : proteintxt = ">" + protein.code1
            if format == 'sp' : proteintxt = ">" + "sp|" + protein.code1 + "|" + protein.code2 + " " + protein.description
            file.write(proteintxt)
            file.write("\n")

            if len(protein.sequence) <= chunksize : 
                file.write(protein.sequence)
                continue
            
            chunk     = protein.sequence[:chunksize]
            remaining = protein.sequence[chunksize:]
            file.write(chunk)
            file.write("\n")
            while True :
                if len(remaining) <= chunksize :
                    file.write(remaining)
                    file.write("\n")
                    break
                chunk = remaining[:chunksize]
                remaining = remaining[chunksize:]
                file.write(chunk)
                file.write("\n")


    def get_proteins_containing_peptide(self, sequence) :
        """It relates the peptide sequence with proteins in the protein dictionary"""
        proteins = []
        
        #Remove any modification from the given sequence
        sequence_to_search = sequence.upper()
        while sequence_to_search.find('[') > -1 :
            mod_left = sequence_to_search.find('[')
            mod_right = sequence_to_search.find(']')
            if mod_right < len(sequence_to_search) :
                sequence_to_search = sequence_to_search[:mod_left-1] + sequence_to_search[mod_right+1:]
            else : sequence_to_search = sequence_to_search[:mod_left-1]
        
        #Search through the protein dictionary
        for protAccession,protein in self.proteinDictionary.iteritems() :
            protSequence = protein.sequence
            if protSequence.find(sequence_to_search) > -1 : proteins.append(protein)
        
        return proteins
    
    
    def readFasta(self, fastaFileName, decoytag = 'DECOY_'):
        
        #check wether fastaFileName exists or not
        if not os.path.exists(fastaFileName):
            print("The file: %s does not exist!" % fastaFileName)
            sys.exit(2)
        
        is_chunksize_set = False
        counter = 0
        with open(fastaFileName,"r") as file:    
            proteinToSend = False
            code1 = ''
            code2 = ''
            modres = ''
            ncbi_tax_id = ''
            description = ''
            sequence = ''
            
            for line in file:
                if len(line)>0:
                    if line[0] =='>' :
                        #send the last protein stored to the database
                        if proteinToSend==True:
                            #~ print (code1, code2, modres, ncbi_tax_id, description) 
                            ##                            if counter ==1 : 
                            ##                                print code1, code2, modres, ncbi_tax_id, description, fastaFileName
                            ##                                print ""
                            ##                                print sequence
                            
                            #remove the absurd astherisk at the end of the sequence that some databases include (sigh)
                            if sequence[-1:] == '*' : sequence = sequence[:-1]
                            
                            proteinDetails = { 'code1' : code1, 'code2' : code2 , 'modres' : modres, 'ncbi_tax_id' : ncbi_tax_id, 'description' : description, 'sequence' : sequence }
                            protein = Protein(proteinDetails)
                            self.proteinDictionary[code1] = protein
                            
                            #Remove info so that it does not accumulate for the next protein
                            code1 = ''
                            code2 = ''
                            modres = ''
                            ncbi_tax_id = ''
                            description = ''
                            sequence = ''
                            
                        
                        #set up to send a new protein to the database
                        dectag = '>' + decoytag
                        if line[:len(dectag)] == dectag : proteinToSend = False
                        if line[:len(dectag)] != dectag : 
                            
                            # Try the first annotation structure
                            firstannotation = False
                            counter += 1
                            proteinToSend = True
                            sequence = ''
                            values = line.split('\\')
                            #if counter == 1: print values
                            for value in values:
                                #'clean' values
                                value = value.strip()
                                if value[-2:] == '\\n': value = value[:-2]
                                
                                #fetch values
                                if value[:1]    == '>'           : code1 = value[1:]
                                elif value[:3]  == 'ID='         : 
                                    code2 = value[3:]
                                    firstannotation = True
                                elif value[:7]  == 'MODRES='     : modres = value[7:]
                                elif value[:10] == 'NCBITAXID='  : ncbi_tax_id = value[10:] 
                                elif value[:3]  == 'DE='         : description = value[3:]
                            
                            #If the first annotation was not fetched, then try the second one.
                            if not firstannotation :
                                values = line.split('|')
                                if len(values) > 2 :
                                    code1 = values[1].strip()
                                    
                                    namesplit = values[2].split(' ')
                                    
                                    code2 = namesplit[0]
                                    
                                    description = ''
                                    for val in namesplit :
                                        if '=' not in val[:] and code2 not in val : 
                                            description += val
                                            description += ' '
                                else : # Third annotation (sigh!)
                                    values = values[0].split(' ')
                                    if len(values) < 2 : continue
                                    if values[0][0] == '>' :
                                        code1 = values[0][1:].strip()
                                        for val in values :
                                            if 'SGDID:' in val :
                                                code2 = val.split(':')[1].strip()
                                                if code2[-1:] == ',' : code2 = code2[:-1]
                                    #So far it does not parse the description, as it is quite complicated to parse... OMG!! Who did design this "nomenclature" ????
                    
                    ##                                First annotation structure
                    ##                            >DECOY_P32527 
                    ##                            \ID=ZUO1_YEAST 
                    ##                            \MODRES=(50|Phosphoserine.)(63|Phosphothreonine.)(114|Phosphoserine.)(426|Phosphoserine.) 
                    ##                            \NCBITAXID=4932 
                    ##                            \DE=Zuotin (J protein ZUO1)
                    
                    ##                                Second annotation structure
                    ##                            >sp|P65637|34KD_MYCTU 34 kDa antigenic protein homolog OS=Mycobacterium tuberculosis GN=Rv0954
                    
                    ##                            Third annotation structure
                    ##                            >YAL001C TFC3 SGDID:S000000001, Chr I from 151006-147594,151166-151097, Genome Release 64-1-1, reverse complement, Verified ORF, "Largest of six subunits of the RNA polymerase III transcription initiation factor complex (TFIIIC); part of the TauB domain of TFIIIC that binds DNA at the BoxB promoter sites of tRNA and similar genes; cooperates with Tfc6p in DNA binding
                    
                    else:
                        #fetch the sequence to the current protein
                        sequence += line[:].rstrip()
                        if not is_chunksize_set : 
                            self.chunksize = len(line[:-1])
                            is_chunksize_set = True
        
            if proteinToSend == True :  # Send the very last protein of the file
                if sequence[-1:] == '*' : sequence = sequence[:].rstrip()
                
                proteinDetails = { 'code1' : code1, 'code2' : code2 , 'modres' : modres, 'ncbi_tax_id' : ncbi_tax_id, 'description' : description, 'sequence' : sequence }
                protein = Protein(proteinDetails)
                self.proteinDictionary[code1] = protein
            
        
        if file.closed: print("file is closed")
        
        print("Number of target proteins: %s" % counter)


def usage() :
    print("proteinDB.py")
    print("-" * 16)
    print("Given a peptide list, it gives back the protein(s) for each peptide.")
    print("Usage: ")
    print("python proteinDB.py -f fastafile.fasta peptidelist.txt")

def readPeptideListCSV(filename):
    
    print("Reading target peptide list csv file %s" % filename)
    reader = csv.reader(open(filename), dialect ='excel-tab')
    
    v= [ row[0] for row in reader if row ]
    
    print("...done.")
    
    return v

def removeModifications(peptides) :
    peptide_output = []
    for peptide in peptides :
        while peptide.find('[') > -1 :
            left_bracket  = peptide.find('[')
            right_bracket = peptide.find(']')
            peptide_tmp = peptide[:left_bracket] + peptide[right_bracket+1:]
            peptide = peptide_tmp
        peptide_output.append(peptide)
            
    return peptide_output
    
def writecsv(headers,matrix,filename) :

    try :
        writer = csv.writer(open(filename,'w'), dialect='excel-tab')
    except :
        print("something went wrong while writing this file : " , filename )
        sys.exit(1)
    
    writer.writerow(headers)
    
    for row in matrix :
        writer.writerow(row)

def main(argv) :
    
    fastaFile  = ''
    outputFile = ''
    
    #Get options
    try:
        opts, args = getopt.getopt(argv, "f:ho:",["fasta=","help","output="])
    
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    argsUsed = 0
    
    for opt,arg in opts:
        #print opt, arg
        if opt in ("-f","--fasta") :
            argsUsed +=2
            fastaFile = arg
        if opt in ("-o","--output") :
            argsUsed +=2
            outputFile = arg 
        elif opt in ("-h", "--help") :
            usage()
            sys.exit()
        
    
    if len(fastaFile) == 0 : 
        print("You must input a fasta file!")
        usage()
        sys.exit()

    peptidesfile = argv[argsUsed] 
    
    #Read fasta file
    protDB = ProteinDB()
    protDB.readFasta(fastaFile)
    
    #Read peptide list
    peptides = readPeptideListCSV(peptidesfile)
    peptides = removeModifications(peptides)
    
    headers = ['peptide','num_references','proteins']
    
    if len(outputFile) == 0 :
        for peptide in peptides :
            print(peptide, len(protDB.findSequenceInProteins(peptide)) , protDB.findSequenceInProteins(peptide))
        sys.exit()
    
    
    
    outputmatrix = []

    for peptide in peptides :
        
        row = []
        row.append(peptide)
        proteins = protDB.findSequenceInProteins(peptide)
        proteins_txt = ''
        row.append(len(proteins))
        for pr in proteins : 
            proteins_txt += pr
            proteins_txt += ','
        
        if proteins_txt[-1:] == ',' : proteins_txt = proteins_txt[:-1]
        
        row.append(proteins_txt)
        outputmatrix.append(row)
        
    writecsv(headers, outputmatrix, outputFile)
    
        
if __name__ == "__main__":    
    main(sys.argv[1:])




