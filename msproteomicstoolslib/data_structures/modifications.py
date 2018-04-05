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
try:
    from elements       import Formulas        
    from aminoacides    import Aminoacides
    from peptide        import Peptide
except ImportError:
    from .elements       import Formulas        
    from .aminoacides    import Aminoacides
    from .peptide        import Peptide

import csv
import ast
import re
import pkg_resources
import sys

class Modifications:
    """
    A collection of modifications
    """
    
    def __init__(self, default_mod_file=None):
        self.list= []
        self.mods_TPPcode = {}   # a more confortable way to store the modifications, if you want just to declare them elsewhere in a sequence
        self.mods_unimods = {}
        self._initModifications(default_mod_file)
    
    def _initModifications(self, default_mod_file):
        if not default_mod_file:
            default_mod_file = pkg_resources.resource_filename(__name__, "modifications_default.tsv")
        self.readModificationsFile(default_mod_file)

    def appendModification(self, modification) :
        self.list.append(modification)
        self.mods_TPPcode[modification.TPP_Mod] = modification
        self.mods_unimods [modification.unimodAccession] = modification
        
    def is_bool(self,expression) :
        return expression.lower() in ['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh', 'por supuesto', 'of course']
    
    def printModifications(self) :
        
        for mymod in self.list:
            print(["%s : %s" % (prop,val) for prop,val in vars(mymod).items() ])
        
    def readModificationsFile(self, modificationsfile):
        '''It reads a tsv file with additional modifications. Modifications will be appended to the default modifications 
        of this class.
        Tsv file headers & an example:
        modified-AA    TPP-nomenclature    Unimod-Accession    ProteinPilot-nomenclature    is_a_labeling    composition-dictionary
        S    S[167]    21    [Pho]    False    {'H' : 1,'O' : 3, 'P' : 1}
        '''

        # sys.stderr.write("reading modifications file: %s \n" % modificationsfile)
        reader = csv.reader(open(modificationsfile,'r'),dialect="excel-tab")
        headers = ['modified-AA', 'TPP-nomenclature',   'Unimod-Accession',  'ProteinPilot-nomenclature', 'is_a_labeling',
                   'composition-dictionary']
         
        header_found = False    
         
        for row in reader:
            if row[0] in headers : 
                header_found = True
                header_d = dict([(l, i) for i, l in enumerate(row)])
                continue
           
            if not header_found:
                continue
        
            mod = Modification(row[header_d['modified-AA']], row[header_d['TPP-nomenclature']], int(row[header_d['Unimod-Accession']]), 
                               row[header_d['ProteinPilot-nomenclature']], self.is_bool(row[header_d['is_a_labeling']]), 
                               ast.literal_eval(row[header_d['composition-dictionary']]) )
                        
            self.appendModification(mod)


    def translateModificationsFromSequence(self, sequence, code, aaLib = None) :
        '''Returns a Peptide object, given a sequence with modifications in any of the available codes.
        The code (TPP, Unimod,...) to be translated must be given.'''
        
        if code not in Modification.codes : 
            #Throw an Exception
            print("The following nomenclature (code) is not recognized : ", code)
            sys.exit(5)
        
        aminoacides_with_mods = []
        terminal_mods = []
        if code == 'unimod' :
            aminoacides_with_mods     = re.findall('([A-Z]\([^\)]*\)|[A-Z])', sequence )
        else :
            aminoacides_with_mods     = re.findall('([A-Z]\[[^\]]*\]|[A-Z])', sequence )
            #Warning : this will only work for TPP unfortunately. It is though the most common operation we'll do
            if code == 'TPP' : terminal_mods = re.findall('([a-z]\[[^\]]*\]|[a-z])', sequence) 
                    
        #mods_peptide is a dictionary wich uses the position of the modification as key, and a Modification object as value:
        #example : GGGGMoxDDCDK  -> mods_peptide = { 5 : Modification1 , 8 : Modification2 }
        mods_peptide = {}
        sequence_no_mods = ''
        for i,aa in enumerate(aminoacides_with_mods) :
            #clean the sequence (basically, take the first letter of the aminoacides_with_mods)
            if len(aa) > 0 : sequence_no_mods += aa[0]
            if len(aa) > 1 : #search for the modification
                modification_found = False
                for modif in self.list :
                    if aa == modif.getcode(code) :
                        #add the modification to the mods_peptide dictionary
                        mods_peptide[i+1] = modif
                        modification_found = True
                if not modification_found :
                    #Throw an Exception
                    print("This modification has not been recognized : " , aa)
                    print("Found in the following sequence : " , sequence)
                    print("The code used to interpret it was : " , code)
                    sys.exit(4) 
        
        for i, mod in enumerate(terminal_mods) :
            modification_found = False
            for modif in self.list :
                if mod == modif.getcode(code) :
                    modification_found = True
                    if modif.is_Nterminal : mods_peptide[0] = modif
                    if modif.is_Cterminal : mods_peptide[len(sequence_no_mods)+1] = modif
            if not modification_found :
                #Throw an Exception
                print("This modification has not been recognized : " , mod)
                print("Found in the following sequence : " , sequence)
                print("The code used to interpret it was : " , code)
                sys.exit(4) 
        
        return Peptide(sequence_no_mods, mods_peptide, aminoacidLib = aaLib)
        
class Modification:
    """
    A modification on an Aminoacid
    """
    
    #: Available modification formats
    codes = ['TPP', 'unimod', 'ProteinPilot']
    
    def __init__(self, aminoacid, tpp_Mod, unimodAccession, peakViewAccession, is_labeling, composition):
        #self.aminoacid = Aminoacid()
        self.aminoacid            = aminoacid
        self.is_Nterminal = False
        self.is_Cterminal = False
        if self.aminoacid == 'N-term' : self.is_Nterminal = True
        if self.aminoacid == 'C-term' : self.is_Cterminal = True
        self.unimodAccession     = unimodAccession
        self.peakviewAccession     = peakViewAccession
        self.TPP_Mod             = tpp_Mod
        self.composition         = composition
        self.is_labeling         = is_labeling
        #if deltaMass        : self.deltamass = deltaMass
        self.deltamass             = Formulas.mass(composition) #self._getMass()
        self.id                 = aminoacid + str(int(round(self.deltamass,0)))

    def getcode(self, code):
        if code not in Modification.codes :
            print("Can't process the requested modification code : " , code)
            print("Available codes are: " , Modification.codes)
            sys.exit(5)
        
        if code == 'TPP' :             return self.TPP_Mod
        if code == 'unimod' :        return "%s(UniMod:%s)" % (self.aminoacid, self.unimodAccession)
        if code == 'ProteinPilot' :    return "%s%s" % (self.aminoacid, self.peakviewAccession)
        

def test(args = []):
    mods = Modifications()

    if len(args) > 0 :
        for arg in args :
            mods.readModificationsFile(arg)

    for mymod in mods.list:
        print(["%s : %s" % (prop,val) for prop,val in vars(mymod).items() ])

        
    #Translate some peptide sequences
    sequences =['PEPTIMEK' , 'PEPTIM[147]EK', 'n[43]PEPTIMEK']
    peptides = []
    for sequence in sequences :
        peptide = mods.translateModificationsFromSequence(sequence, "TPP")
        print("peptide sequence : " , peptide.sequence)
        print("peptide modifications :" )
        for mod in peptide.modifications.values() :
            print(mod.id , mod.deltamass)
        print("peptide mass : " , peptide.mass)

        

if __name__ == "__main__":
    import sys
    test(sys.argv[1:])
    sys.exit(2)
    
