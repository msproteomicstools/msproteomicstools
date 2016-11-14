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
    from aminoacides import Aminoacides
    from elements import Elements
    from elements import Formulas
except ImportError:
    from .aminoacides import Aminoacides
    from .elements import Elements
    from .elements import Formulas

#from spectrum import Spectrum
from types import *
import random
import sys
import itertools

class Peptide:
    def __init__(self, sequence,  modifications={}, protein="", aminoacidLib = None): 
        self.sequence = sequence
        #self.modifications is a dictionary wich uses the position of the modification as key, and a Modification object as value:
        #example : GGGGMoxDDCDK  -> self.modifications = { 5 : ModificationObject1 , 8 : ModificationObject2 }
        self.aaList = aminoacidLib
        if self.aaList == None : self.aaList = Aminoacides()
        self.modifications = {}
        if len(modifications) > 0 : self.modifications = modifications
        self.composition = self._getComposition()
        self.sequenceMods = self.getSequenceWithMods("TPP")
        self.mass = self._getMassFromSequence()
        self.spectra = []
        self.proteins = protein
        self.labelings = ['N15', '15N', 'AQUA_KR', 'SILAC_K6R10', 'no_labeling', 'SILAC_K8R10', 'SILAC_K8R6']
        self.iontypes = ['a', 'b', 'c', 'x', 'y', 'z']
        self.lossConstraints = { -97.976896 : ['S','T','Y'], -79.966331 : ['S','T','Y'] }

    def _it(self,l) :
        try : iter(l)
        except TypeError : iterable = [l]
        else : iterable = l
        return iterable

    def all_ions(self, ionseries = None, frg_z_list = [1,2], fragmentlossgains = [0,] , mass_limits = None, label = '' ):
        '''
        Returns all the fragment ions of the peptide in a tuple of two objects: (annotated, ionmasses_only) 
        annotated is a list of tuples as : (ion_type, ion_number, ion_charge, lossgain, fragment_mz)
        ionmasses_only is a list of fragment masses.
        When ionseries is not provided, all existing ion series (see: Peptide.iontypes) will be calculated.
        When frg_z_list is not provided, fragment ion charge states +1 and +2 will be used.
        '''
        annotated = []
        ionmasses_only = []
        
        if not ionseries : ionseries = self.iontypes
        
        for ion_type in ionseries :
            for ion_number in range(1,len(self.sequence)) :
                for ion_charge in frg_z_list :
                    for lossgain in fragmentlossgains :
                        frg_seq = self.fragmentSequence(ion_type, ion_number)
                        lossgainOK = True
                        for mass, constrain in self.lossConstraints.items() :
                            if abs(mass-lossgain) < 0.02 :
                                lossgainOK = False 
                                for aa in constrain : 
                                    if aa in frg_seq : 
                                        lossgainOK = True
                                        break
                        if not lossgainOK : continue
                                
                        frg_mz = self.getMZfragment(ion_type, ion_number, ion_charge, label='', fragmentlossgain = lossgain)
                        if mass_limits : 
                            if frg_mz < mass_limits[0] : continue
                            if frg_mz > mass_limits[1] : continue
                        annotated.append( ( ion_type, ion_number, ion_charge, lossgain, frg_mz ) )
                        ionmasses_only.append(frg_mz)
        
        return ( annotated, ionmasses_only )

        
    def cal_UIS(self, otherPeptidesList, UISorder = 2,  ionseries = None, fragmentlossgains = [0,], precision = 1e-8, frg_z_list = [1,2], mass_limits = None) :
        '''
        It calculates the UIS for a given peptide referred to a given list of other peptides.
        It returns a tuple of two objects all_UIS, and all_UIS_annotated. all_UIS contains only a mass list.
        '''
        if not ionseries : ionseries = self.iontypes
        
        all_UIS = [] #   all_UIS = [ [UIS1] , [UIS2] ,[UIS3]  ]
        all_UIS_annotated = {} # all_UIS_annotated = { UIS1 : [annotation(s)] , ... }
        print (fragmentlossgains)
        selfPep_annotated , selfPep_masses = self.all_ions( ionseries = ionseries, frg_z_list = frg_z_list, 
                                                        fragmentlossgains = fragmentlossgains , mass_limits = mass_limits, 
                                                        label = '' )
        if UISorder < 0 : 
            all_UIS = [selfPep_masses]
            all_UIS_annotated = { tuple(selfPep_masses) : selfPep_annotated }
            return (all_UIS, all_UIS_annotated)
        
        ions_others = []  # ions_others = [ [ions_pep1] , [ions_pep2], ]
        #Create a pool of the ions of other peptides         
        for othPep in otherPeptidesList :
            _ , othIons = othPep.all_ions(ionseries = ionseries, frg_z_list = frg_z_list, fragmentlossgains = fragmentlossgains , 
                            mass_limits = mass_limits, label = '')
            ions_others.append(othIons)
        
        #Remove from the selfPep_annotated and selfPep_masses the ions that are common to ALL the other peptides 
        #(since they are not going to be useful at all!)
        discard = []
        for ion in selfPep_masses : 
            in_all = True
            for pep in ions_others : 
                if not set(self._it(ion)).issubset(pep) : in_all = False
            if in_all : discard.append(ion)
        selfPep_masses = [mass for mass in selfPep_masses if mass not in discard ]
        discard = []
        #Starting for the lowest order of UIS, so that we can remove combinations containing a unique sub-combination
        for order in range(1,UISorder+1) :
            tentative_UIS = list(itertools.combinations(selfPep_masses, order))
            #Remove from the tentative_UIS the sub-combinations that are already UIS
            for uis in all_UIS :
                for t_uis in tentative_UIS :
                    if set(uis).issubset(t_uis) :  discard.append(t_uis)
            #print "discard : " , discard
            tentative_UIS_tmp = [ tuis for tuis in tentative_UIS if not set(tuis).issubset(discard)]
            tentative_UIS = tentative_UIS_tmp
            discard = []
            for i, t_uis in enumerate(tentative_UIS) :
                #Check whether the tentative UIS actually defines uniquely the peptide
                for othPep in ions_others :
                    if set(t_uis).issubset(othPep) :
                        #print "removing : " , t_uis , tentative_UIS[i] 
                        discard.append(t_uis)
            tentative_UIS_tmp = [ tuis for tuis in tentative_UIS if tuis not in discard]
            tentative_UIS = tentative_UIS_tmp
            
            for t_uis in tentative_UIS :
                all_UIS.append(t_uis)
                #Annotate the UIS, and append it to the all_UIS_annotated
                uis_annotations = []
                for mass in t_uis : 
                    for i, annot in enumerate(selfPep_annotated) :
                        if mass == annot[4] : uis_annotations.append(annot)
                        #break
                all_UIS_annotated[t_uis] = uis_annotations

        return ( all_UIS , all_UIS_annotated )

    def comparePeptideFragments(self, otherPeptidesList, ionseries = None, fragmentlossgains = [0,] , precision = 1e-8, frg_z_list = [1,2]) :
        '''
        This returns a tuple of lists: (CommonFragments, differentialFragments). The differentialFragmentMasses are the masses
        of the __self__ peptide are not shared with any of the peptides listed in the otherPeptidesList. otherPeptidesList must be a list of
        Peptide objects. The fragments are reported as a tuple : (ionserie,ion_number,ion_charge,frqgmentlossgain,mass)
        '''
        if not ionseries : ionseries = self.iontypes
        
        ions_selfPep = []
        ions_others  = [] 
        
        for ion_type in ionseries :
            for ion_number in range(1,len(self.sequence)) :
                for ion_charge in frg_z_list :
                    for lossgain in fragmentlossgains :
                        ions_selfPep.append( ( ion_type, ion_number, ion_charge, lossgain, self.getMZfragment(ion_type, ion_number, ion_charge, label='', fragmentlossgain = lossgain) ) ) 
        
        for othPep in otherPeptidesList :
            for ion_type in ionseries :
                for ion_number in range(1,len(othPep.sequence)) :
                    for ion_charge in frg_z_list :
                        for lossgain in fragmentlossgains :
                            ions_others.append(( ion_type, ion_number, ion_charge, lossgain, othPep.getMZfragment(ion_type, ion_number, ion_charge, label='', fragmentlossgain=lossgain) ) )
        
        shared_ions = []
        unmatched_ions = [] 
        for (ion_type, ion_number, ion_charge, fragmentlossgain, mass1) in ions_selfPep :
            is_mass_match = False
            for (_,_,_,_,mass2) in ions_others :
                if abs(mass1-mass2) <= precision : 
                    shared_ions.append((ion_type, ion_number, ion_charge, fragmentlossgain, mass1))
                    is_mass_match = True
            if not is_mass_match : unmatched_ions.append((ion_type, ion_number, ion_charge, fragmentlossgain, mass1))
        
        return (shared_ions, unmatched_ions)
    
    def calIsoforms(self, switchingModification, modLibrary) :
        '''This returns the full list of peptide species of the same peptide family (isobaric, same composition, different 
        modification site. The list is given as a list of Peptide objects. switchingModification must be given as a Modification object.
        '''
        peptide_family = []
        #Count the number of switching modifications 
        num_of_switchMods = 0
        fixed_modifications = {}
        for key, mod in self.modifications.items() :
            if mod.unimodAccession == switchingModification.unimodAccession : num_of_switchMods += 1
            else : fixed_modifications[key] = mod
            
        #Store the position of the available modification sites of the peptide
        modSites = []
        for site,aa in enumerate(self.sequence[:]) :
            for mod in modLibrary.list :
                #print mod.unimodAccession , switchingModification.unimodAccession 
                if mod.unimodAccession == switchingModification.unimodAccession :
                    if mod.aminoacid == aa : modSites.append(site)  
        
        for subset in itertools.combinations(modSites, num_of_switchMods) :
            newmods = dict([ (key, mod) for key, mod in fixed_modifications.items() ])  #Keep the fixed modifications!
            for idx in subset : newmods[idx+1] = switchingModification
            isoform = Peptide(self.sequence, newmods, aminoacidLib = self.aaList)
            peptide_family.append(isoform)
        
        return peptide_family

    def getSequenceWithMods(self, code) :
        seqMods = ''
        
        for i, aa in enumerate([''] + list(self.sequence) + [''], start=-1):
            if i + 1 in self.modifications:
                codestr = self.modifications[i + 1].getcode(code)
                code_sans_aa = codestr[len("*-term"):] if codestr.startswith(("N-term", "C-term")) else \
					codestr if codestr == "n[43]" else \
						codestr[1:]
                seqMods += aa + code_sans_aa
            else : seqMods += aa
            
        return seqMods


    def _getMassFromSequence(self):
        #self.aaList = Aminoacides()

        mass = 0
        for aa in self.sequence[:]:
            for ac in self.aaList.list:
                if aa == ac.code:
                    mass += ac.deltaMass
        
        #Adding an H2O to the mass, as it is calculated by using delta masses
        mass += Formulas.mass(Formulas.H2O) #(1.007825032 *2 + 15.99491462) 
        
        #Adding modifications deltamass
        for aaPosMod in self.modifications.keys() :
            mass += self.modifications[aaPosMod].deltamass

        return mass
    
        
    def getDeltaMassFromSequence(self, sequence):
        
        #self.aaList = Aminoacides()
        
        mass = 0
        for aa in sequence[:]:
            for ac in self.aaList.list:
                if aa == ac.code:
                    mass += ac.deltaMass
       
        return mass

    def pseudoreverse(self, sequence='None'):
        
        if sequence == 'None':        
            revseq = self.sequence[::-1][1:] + self.sequence[-1:]
            
            self.sequence = revseq
        else:
            revseq = sequence[::-1][1:] + sequence[-1:]
            
            return revseq

    def shuffle_sequence(self):

        list = []
        for aa in self.sequence[:-1]:
            list.append(aa)
            
        random.shuffle(list)

        sequence_shuffled = ''

        for aa in list:
            sequence_shuffled += aa
            
        sequence_shuffled += self.sequence[len(self.sequence) - 1]

        self.sequence = sequence_shuffled
        
    def get_decoy_Q3(self, frg_serie, frg_nr, frg_z, blackList=[], max_tries=1000):
        
        old_sequence = self.sequence

        #If b- ion and frg-nr is the second to last, it wont be impossible to get a different mass -> we use other frg_nr
        if frg_serie in ['b'] and int(frg_nr) > len(self.sequence) - 1:
            frg_nr = str (int(frg_nr) - 1)
        
        try_ = 1
        found = False
        notvalid = False
        q3_decoy = 0.0
        while (try_ < max_tries and found == False):
            notvalid = False
            self.shuffle_sequence()
            q3_decoy = self.getMZfragment(frg_serie, frg_nr, frg_z)
            for q3 in blackList:
                #print abs(q3-q3_decoy)
                if abs(q3 - q3_decoy) < 0.001:
                    #print q3
                    notvalid = True
            if notvalid == False: found = True

            try_ += 1
            
            
        self.sequence = old_sequence
        
        return q3_decoy
    
    def _getComposition(self):
        #aaList = Aminoacides()
        composition = {}
        
        for aa in self.sequence[:]:
            for ac in self.aaList.list:
                if aa == ac.code:
                    for elem, num in ac.composition.items():
                        if elem not in composition:
                            composition[elem] = num
                        else:
                            composition[elem] += num
                        
        for aaModified, modification in self.modifications.items() :
            for elem,numAtoms in modification.composition.items() :
                if elem not in composition : composition[elem] = numAtoms
                else : composition[elem] += numAtoms
        
        return composition
    
    def _getCompositionSeq(self, sequence, modifications = {}):
        #aaList = Aminoacides()
        composition = {}
        
        for aa in sequence[:]:
            for ac in self.aaList.list:
                if aa == ac.code:
                    for elem, num in ac.composition.items():
                        if elem not in composition:
                            composition[elem] = num
                        else:
                            composition[elem] += num

        for modification in modifications :
            for elem,numAtoms in modification.composition.items() :
                if elem not in composition : composition[elem] = numAtoms
                else : composition[elem] += numAtoms
        
        return composition
        
    
    def _getAminoacidList(self, fullList=False):
        
        aminoacides = self.aaList
        aaList = {}
        
        if fullList:            
            for ac in aminoacides.list:
                aaList[ac.code] = 0

        for aa in self.sequence[:]:
            for ac in aminoacides.list:
                if aa == ac.code:
                    if aa not in aaList:
                        aaList[aa] = 1
                    else:
                        aaList[aa] += 1
            
        return aaList


    def getMZ(self, charge, label=''):
        
        #Check label
        if len(label) > 0 and label not in self.labelings : 
            print("Coding error: this labeling is not reported!! %s" % label)
            sys.exit(2)
            
        
        pepmass = self._getMassFromSequence()
        z = charge
        H = 1.007825032
        massC = 12.0000000
        massProton = 1.0072765
        massN = 14.00307401
        massC13 = 13.003355
        massN15 = 15.0001089
        massShiftN15 = massN15 - massN
        massShiftC13 = massC13 - massC
        
        
        
        # N15 labeling: adds up the shift mass
        shift = 0.0
        if label == '15N' or label == 'N15' :
            composition = self._getComposition()
            if 'N' in composition:
                shift = (float(composition['N']) * massShiftN15) 
                
        #AQUA_KR labeling: adds up the shift mass
        if label == 'AQUA_KR':
            #Update (23.nov.2010): AQUA_KR only modifies Lys and Arg in the C-terminal position.  
            if 'K' in self.sequence[-1:] : shift += 6 * massShiftC13 + 2 * massShiftN15
            if 'R' in self.sequence[-1:] : shift += 6 * massShiftC13 + 4 * massShiftN15

        if label == 'SILAC_K6R10' :
            #All K and R amino acides are shifted 6 Da and 10 Da.
            numK = 0
            numR = 0
            for aa in self.sequence[:] :
                if aa == 'K' : numK += 1
                if aa == 'R' : numR += 1
            shift += numK * (6 * massShiftC13) 
            shift += numR * (6 * massShiftC13 + 4 * massShiftN15)

        if label == 'SILAC_K8R10' :
            #All K and R amino acides are shifted 8 Da and 10 Da.
            numK = 0
            numR = 0
            for aa in self.sequence[:] :
                if aa == 'K' : numK += 1
                if aa == 'R' : numR += 1
            shift += numK * (6 * massShiftC13 + 2 * massShiftN15) 
            shift += numR * (6 * massShiftC13 + 4 * massShiftN15)
            
        if label == 'SILAC_K8R6' :
            #All K and R amino acides are shifted 8 and 6 Da.
            numK = 0
            numR = 0
            for aa in self.sequence[:] :
                if aa == 'K' : numK += 1
                if aa == 'R' : numR += 1
            shift += numK * (6 * massShiftC13 + 2 * massShiftN15)
            shift += numR * (6 * massShiftC13)

 
        mz = (pepmass + shift + massProton * z) / z
       
        return mz

    def fragmentSequence(self, ion_type, frg_number) :
        if ion_type == 'p' :
            return self.sequence
        elif ion_type == 'a' :
            return self.sequence[:frg_number] 
        elif ion_type == 'b' :
            return self.sequence[:frg_number] 
        elif ion_type == 'c' :
            return self.sequence[:frg_number] 
        elif ion_type == 'x' :
            return self.sequence[-frg_number:]
        elif ion_type == 'y' :
            return self.sequence[-frg_number:]
        elif ion_type == 'z' :
            return self.sequence[-frg_number:]

        return self.sequence
    
    def getMZfragment(self, ion_type, ion_number, ion_charge, label='', fragmentlossgain=0.0):
        #Check label
        if len(label) > 0 and label not in self.labelings : 
            print("Coding error: this labeling is not reported!! %s" % label)
            sys.exit(2)
        
        #Check ion type
        #if ion_type not in self.iontypes :
            #print "This ion type can not be processed :" , ion_type
            #print "This is the list of ion types considered in this library : " , self.iontypes
            #sys.exit(2)
        
        massH = 1.007825032
        massO = 15.99491462
        massN = 14.00307401
        massC = 12.0000000
        massProton = 1.0072765
        massNH3 = Formulas.mass(Formulas.NH3) #massN + 3 * massH
        massH2O = Formulas.mass(Formulas.H2O) #(massH *2 + massO) 
        massCO2     = Formulas.mass(Formulas.CO2) #(massC + 2 * massO)
        protonMass = 1.007825032

        massN15 = 15.0001089
        massShiftN15 = massN15 - massN
        massC13 = 13.003355
        massShiftC13 = massC13 - massC
       
        frg_number = int(ion_number)
        frg_charge = int(ion_charge)    
        
        mzfragment = 0
        frg_seq = ''

        #precursors
        if ion_type == 'p' :
            mzfragment = self.getMZ(frg_charge, label)
            return mzfragment
        
        #a series
        if ion_type == 'a' :
            frg_seq = self.sequence[:frg_number] 
            mzfragment = self.getDeltaMassFromSequence(frg_seq)
            mzfragment -= (massC + massO)
            #Adding modifications deltamass
            for aaPosMod in self.modifications.keys() :
                if aaPosMod <= frg_number: mzfragment += self.modifications[aaPosMod].deltamass

        #b series
        if ion_type == 'b' :
            frg_seq = self.sequence[:frg_number] 
            mzfragment = self.getDeltaMassFromSequence(frg_seq)
            #Adding modifications deltamass
            for aaPosMod in self.modifications.keys() :
                if aaPosMod <= frg_number: mzfragment += self.modifications[aaPosMod].deltamass
        
        #c series
        if ion_type == 'c' :
            frg_seq = self.sequence[:frg_number] 
            mzfragment = self.getDeltaMassFromSequence(frg_seq)
            mzfragment += massNH3
            #Adding modifications deltamass
            for aaPosMod in self.modifications.keys() :
                if aaPosMod <= frg_number: mzfragment += self.modifications[aaPosMod].deltamass
        
        #x series
        if ion_type == 'x' :
            frg_seq = self.sequence[-frg_number:]
            mzfragment = self.getDeltaMassFromSequence(frg_seq)
            mzfragment += massCO2
            #Adding modifications deltamass
            for aaPosMod in self.modifications.keys() :
                if (len(self.sequence) - frg_number) < aaPosMod : mzfragment += self.modifications[aaPosMod].deltamass


        #y series
        if ion_type == 'y' :
            frg_seq = self.sequence[-frg_number:]
            mzfragment = self.getDeltaMassFromSequence(frg_seq)
            mzfragment += massH2O
            #Adding modifications deltamass
            for aaPosMod in self.modifications.keys() :
                if (len(self.sequence) - frg_number) < aaPosMod : mzfragment += self.modifications[aaPosMod].deltamass

        #z series
        if ion_type == 'z' :
            frg_seq = self.sequence[-frg_number:]
            mzfragment = self.getDeltaMassFromSequence(frg_seq)
            mzfragment += massH2O
            mzfragment -= massNH3
            #Adding modifications deltamass
            for aaPosMod in self.modifications.keys() :
                if (len(self.sequence) - frg_number) < aaPosMod : mzfragment += self.modifications[aaPosMod].deltamass

            
        # 15N labeling: adds up the shift mass
        if label == '15N' or label == 'N15':
            frg_composition = self._getCompositionSeq(frg_seq)
            if 'N' in frg_composition:
                shift = (float(frg_composition['N']) * massShiftN15) 
                mzfragment += shift
                
        #AQUA_KR labeling: adds up the shift mass
        if label == 'AQUA_KR':
            #Update (23.nov.2010): AQUA_KR only modifies Lys and Arg in the C-terminal position. 
            #Update (1.dec.2010) : I am a jerk. I should have taken it only for y ions , and not for the b ones 
            if 'K' in self.sequence[-1:] and ion_type == 'y' : mzfragment += 6 * massShiftC13 + 2 * massShiftN15
            if 'R' in self.sequence[-1:] and ion_type == 'y' : mzfragment += 6 * massShiftC13 + 4 * massShiftN15
        
        if label == 'SILAC_K6R10' :
            #All K and R amino acides are shifted 6 Da and 10 Da.
            numK = 0
            numR = 0
            for aa in frg_seq[:] :
                if aa == 'K' : numK += 1
                if aa == 'R' : numR += 1
            mzfragment += numK * (6 * massShiftC13)
            mzfragment += numR * (6 * massShiftC13 + 4 * massShiftN15)
        
        if label == 'SILAC_K8R10' :
            #All K and R amino acides are shifted 8 Da and 10 Da.
            numK = 0
            numR = 0
            for aa in frg_seq[:] :
                if aa == 'K' : numK += 1
                if aa == 'R' : numR += 1
            mzfragment += numK * (6 * massShiftC13 + 2 * massShiftN15) 
            mzfragment += numR * (6 * massShiftC13 + 4 * massShiftN15)

        if label == 'SILAC_K8R6' :
            #All K and R amino acides are shifted 8 Da and 6 Da.
            numK = 0
            numR = 0
            for aa in frg_seq[:] :
                if aa == 'K' : numK += 1
                if aa == 'R' : numR += 1
            mzfragment += numK * (6 * massShiftC13 + 2 * massShiftN15) 
            mzfragment += numR * (6 * massShiftC13)


        mzfragment += frg_charge * massProton
        mzfragment += fragmentlossgain
        mzfragment /= frg_charge
        
        return mzfragment
    
    
    
    def addSpectrum(self, spectrum):
        '''Deprecated definition'''
        if isinstance(spectrum, Spectrum):
            self.spectra.append(spectrum)
        #else: #catch error


def test():

    try:
        from modifications import Modifications
    except ImportError:
        from .modifications import Modifications


    mods = Modifications()
    mypep2  = Peptide('LMGPTSVVMGR', modifications={ 9: mods.mods_TPPcode['M[147]']}) #M[147]
    mypep = Peptide('LMGPTSVVMGR')
    phospho    = mods.mods_unimods[21]
    oxi        = mods.mods_unimods[35]
    
    isoform4_1  = Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : oxi , 5 : phospho })
    isoform4_2  = Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : oxi , 12 : phospho })
    isoform4_3  = Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : oxi , 13 : phospho })
    isoform4_4  = Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : oxi , 14 : phospho })
    
    isof_target = isoform4_2
    theOtherIsoforms = [ isoform4_1 , isoform4_2, isoform4_3 , isoform4_4 ]
    theOtherIsoforms.remove(isof_target)
    
    for isof in theOtherIsoforms : 
        print(isof.getSequenceWithMods('unimod'))
        annot , ion_masses =  isof.all_ions(ionseries = ['y',], frg_z_list = [1,], fragmentlossgains = [0,] , mass_limits = [300,1500], label = '' )
        
    print(isof_target.getSequenceWithMods('unimod'))
    annot , ion_masses = isof_target.all_ions(ionseries = ['y',], frg_z_list = [1,], fragmentlossgains = [0,] , mass_limits = [300,1500], label = '' )
    
    UISmass , UISannot = isof_target.cal_UIS(theOtherIsoforms, UISorder = 2,  ionseries = ['y',], fragmentlossgains = [0,], precision = 1e-5, frg_z_list = [1,], mass_limits = [300,1500])
    for UIS in UISannot.values() :
        print (UIS)
    
    matched , unmatched = mypep.comparePeptideFragments([mypep2,], ['y','b'], precision = 1e-2)
    print("Matched ions : " , matched)
    print("Unmatched ions : ", unmatched)
    
    my_isoforms = mypep2.calIsoforms(mods.mods_TPPcode['M[147]'], mods)
    print("isoforms of %s" % mypep2.getSequenceWithMods('unimod') )
    for iso in my_isoforms :
        print (iso.getSequenceWithMods('unimod'))
    
    for pep in [mypep, mypep2] : 
            
        print(pep.sequence, pep.mass, pep.modifications)
        print(pep._getComposition())
        print(pep._getAminoacidList(True))
        print(pep.getSequenceWithMods("ProteinPilot"))
        print(pep.getMZ(1))
             
        print('#,' , ", ".join(pep.iontypes))
        masses = []
        for v in range(1, len(pep.sequence) + 1) :
            masses.append(str(v))
            for ionserie in pep.iontypes :
                masses.append(str(pep.getMZfragment(ionserie, v, 1)))
            print(", ".join(masses))
            masses = []
            
            #for property, value in vars(pep).items() :
            #    print property , " : " , value
        
        
        
        '''
        LIGPTSVVMGR 1144.62742895 {9: 15.99491462}
        {'H': 86, 'C': 49, 'S': 1, 'O': 13, 'N': 14}
        {'A': 0, 'C': 0, 'E': 0, 'D': 0, 'G': 2, 'F': 0, 'I': 1, 'H': 0, 'K': 0, 'M': 1, 'L': 1, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': 1, 'T': 1, 'W': 0, 'V': 2, 'Y': 0}
        1145.63470545
        #, a, b, c, x, y, z
        1, 86.096425862, 114.091340482, 131.117889588, 201.098216784, 175.118952228, 158.092403122
        2, 199.180489844, 227.175404464, 244.20195357, 258.11968051, 232.140415954, 215.113866848
        3, 256.20195357, 284.19686819, 301.223417296, 405.155079738, 379.175815182, 362.149266076
        4, 353.254717424, 381.249632044, 398.27618115, 504.223493656, 478.2442291, 461.217679994
        5, 454.302395898, 482.297310518, 499.323859624, 603.291907574, 577.312643018, 560.286093912
        6, 541.334424308, 569.329338928, 586.355888034, 690.323935984, 664.344671428, 647.318122322
        7, 640.402838226, 668.397752846, 685.424301952, 791.371614458, 765.392349902, 748.365800796
        8, 739.471252144, 767.466166764, 784.49271587, 888.424378312, 862.445113756, 845.41856465
        9, 886.506651372, 914.501565992, 931.528115098, 945.445842038, 919.466577482, 902.440028376
        10, 943.528115098, 971.523029718, 988.549578824, 1058.52990602, 1032.55064146, 1015.52409236
        11, 1099.62922614, 1127.62414076, 1144.65068987, 1171.61397, 1145.63470545, 1128.60815634
        LIGPTSVVMGR 1128.63251433 {}
        {'H': 86, 'C': 49, 'S': 1, 'O': 13, 'N': 14}
        {'A': 0, 'C': 0, 'E': 0, 'D': 0, 'G': 2, 'F': 0, 'I': 1, 'H': 0, 'K': 0, 'M': 1, 'L': 1, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': 1, 'T': 1, 'W': 0, 'V': 2, 'Y': 0}
        1129.63979083
        #, a, b, c, x, y, z
        1, 86.096425862, 114.091340482, 131.117889588, 201.098216784, 175.118952228, 158.092403122
        2, 199.180489844, 227.175404464, 244.20195357, 258.11968051, 232.140415954, 215.113866848
        3, 256.20195357, 284.19686819, 301.223417296, 389.160165118, 363.180900562, 346.154351456
        4, 353.254717424, 381.249632044, 398.27618115, 488.228579036, 462.24931448, 445.222765374
        5, 454.302395898, 482.297310518, 499.323859624, 587.296992954, 561.317728398, 544.291179292
        6, 541.334424308, 569.329338928, 586.355888034, 674.329021364, 648.349756808, 631.323207702
        7, 640.402838226, 668.397752846, 685.424301952, 775.376699838, 749.397435282, 732.370886176
        8, 739.471252144, 767.466166764, 784.49271587, 872.429463692, 846.450199136, 829.42365003
        9, 870.511736752, 898.506651372, 915.533200478, 929.450927418, 903.471662862, 886.445113756
        10, 927.533200478, 955.528115098, 972.554664204, 1042.5349914, 1016.55572684, 999.529177738
        11, 1083.63431152, 1111.62922614, 1128.65577525, 1155.61905538, 1129.63979083, 1112.61324172
        '''
   

if __name__ == "__main__":
    test()
    sys.exit(2)
