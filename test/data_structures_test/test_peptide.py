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

from __future__ import print_function
import unittest
import os

import msproteomicstoolslib.data_structures.peptide as peptide
import msproteomicstoolslib.data_structures.modifications as modifications

# from msproteomicstoolslib.data_structures.peptide import *

class TestUnitPeptide(unittest.TestCase):

    def setUp(self):
        self.mods         = modifications.Modifications()
        self.phospho     = self.mods.mods_unimods[21]
        self.oxi        = self.mods.mods_unimods[35]
        self.carbamidomethyl = self.mods.mods_unimods[4]
        self.carbamyl = self.mods.mods_unimods[5]
        
        self.mypep         = peptide.Peptide('LIGPTSVVMGR', modifications={ 9: self.mods.mods_TPPcode['M[147]'] }) #M[147] 
        self.mypep2        = peptide.Peptide('LIGPTSVVMGR')
        self.mypep3        = peptide.Peptide('ELVISLIVESS', modifications = { 5 : self.phospho , 10 : self.phospho })
        self.isoform1      = peptide.Peptide('ELVISLIVESS', modifications = { 5 : self.phospho , 11 : self.phospho })
        self.isoform2      = peptide.Peptide('ELVISLIVESS', modifications = { 10 : self.phospho , 11 : self.phospho })
        self.mypep4        = peptide.Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : self.oxi , 5 : self.phospho })
        self.isoform4_1  = peptide.Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : self.oxi , 5 : self.phospho })
        self.isoform4_2  = peptide.Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : self.oxi , 12 : self.phospho })
        self.isoform4_3  = peptide.Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : self.oxi , 13 : self.phospho })
        self.isoform4_4  = peptide.Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 1 : self.oxi , 14 : self.phospho })
        self.notisoform4 = peptide.Peptide('MHGGTGFAGIDSSSPEVK', modifications = { 14 : self.phospho })
        self.mypep5        = peptide.Peptide('CNTPTYCDLGK', modifications = { 0 : self.carbamyl , 1 : self.carbamidomethyl , 7 : self.carbamidomethyl })
        
    def assertAlmostIn(self, tpl, other, msg, epsilon=1e-10):
        """
        Check whether the given tpl is in the list of tuples "other"
        """
        any_tuple_equal = False
        for other_tpl in other:
            tuple_equal = True
            for a,b in zip(sorted(tpl), sorted(other_tpl)):
                if abs(a-b) > epsilon:
                    tuple_equal = False
            any_tuple_equal = any_tuple_equal | tuple_equal
        self.assertTrue(any_tuple_equal, msg)

    def test_create_peptide(self):
        p = peptide.Peptide('PEPTIDE')

    def test_peptide_unmodified(self):

        # check that the total mass, mz and composition is correct
        pep = peptide.Peptide('LIGPTSVVMGR')
        self.assertAlmostEqual(pep.mass, 1128.63251433)
        self.assertEqual(pep._getComposition(), {'H': 86, 'C': 49, 'S': 1, 'O': 13, 'N': 14})
        self.assertAlmostEqual(pep.getMZ(1), 1129.63979083)
        
        # TODO also check some ion series
        ## print '#,' , ", ".join(pep.iontypes)
        ## masses = []
        ## for v in range(1,len(pep.sequence)+1) :
        ##     masses.append(str(v))
        ##     for ionserie in pep.iontypes :
        ##         masses.append(str(pep.getMZfragment(ionserie,v,1)))
        ##     print ", ".join(masses)
        ##     masses = []
        ##     
        ##     #for property, value in vars(pep).iteritems() :
        ##     #    print property , " : " , value

    def test_peptide_with_modifications(self):

        # check that the total mass, mz and composition is correct
        # pep = peptide.Peptide('LIGPTSVVMGR', modifications={ 9: self.mods.mods_TPPcode['R[166]'] }) #M[147] 
        pep = peptide.Peptide('LIGPTSVVMGR', modifications={ 9: self.mods.mods_TPPcode['M[147]'] }) #M[147] 
        self.assertAlmostEqual(pep.mass, 1144.62742895)
        self.assertEqual(pep._getComposition(), {'H': 86, 'C': 49, 'S': 1, 'O': 14, 'N': 14} )
        self.assertAlmostEqual(pep.getMZ(1), 1145.63470545)
        
        # TODO also check some ion series
        ## print '#,' , ", ".join(pep.iontypes)
        ## masses = []
        ## for v in range(1,len(pep.sequence)+1) :
        ##     masses.append(str(v))
        ##     for ionserie in pep.iontypes :
        ##         masses.append(str(pep.getMZfragment(ionserie,v,1)))
        ##     print ", ".join(masses)
        ##     masses = []
        ##     
        ##     #for property, value in vars(pep).iteritems() :
        ##     #    print property , " : " , value

    def test_peptide_with_modifications_2(self):
          mypep  = peptide.Peptide('LMGPTSVVMGR', modifications={ 9: self.mods.mods_TPPcode['M[147]'] , 2 : self.mods.mods_TPPcode['M[147]']}) #M[147]

    ###################################
    ## All function unit tests      ##
    ###################################

    def test_getSequenceWithMods(self):
        self.assertEqual(self.mypep.getSequenceWithMods("TPP"), 'LIGPTSVVM[147]GR')
        self.assertEqual(self.mypep.getSequenceWithMods("unimod"), 'LIGPTSVVM(UniMod:35)GR')
        self.assertEqual(self.mypep.getSequenceWithMods("ProteinPilot"), 'LIGPTSVVM[Oxi]GR')
        self.assertEqual(self.mypep5.getSequenceWithMods("TPP"), 'n[43]C[160]NTPTYC[160]DLGK')
        self.assertEqual(self.mypep5.getSequenceWithMods("unimod"), '(UniMod:5)C(UniMod:4)NTPTYC(UniMod:4)DLGK')
        self.assertEqual(self.mypep5.getSequenceWithMods("ProteinPilot"), '[CRM]C[CAM]NTPTYC[CAM]DLGK')

    def test_calIsoforms(self):
        self.isoforms = self.mypep3.calIsoforms(self.phospho, self.mods)
        self.isoforms_seqmods = [ isoform.getSequenceWithMods('unimod') for isoform in self.isoforms ]
        self.assertEqual(len(self.isoforms), 3)
        self.assertIn(self.mypep3.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Original peptide not in the isoforms list")
        self.assertIn(self.isoform1.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Isoform1 not in the isoforms list")
        self.assertIn(self.isoform2.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Isoform2 not in the isoforms list")    
        
        self.isoforms = self.mypep4.calIsoforms(self.phospho, self.mods)
        self.isoforms_seqmods = [ isoform.getSequenceWithMods('unimod') for isoform in self.isoforms ]
        self.assertEqual(len(self.isoforms), 4, "Number of estimated isoforms doesn't match!")
        self.assertIn(self.mypep4.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Original peptide not in the isoforms list")
        self.assertIn(self.isoform4_1.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Isoform1 not in the isoforms list")
        self.assertIn(self.isoform4_2.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Isoform2 not in the isoforms list")    
        self.assertIn(self.isoform4_3.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Isoform3 not in the isoforms list")
        self.assertIn(self.isoform4_4.getSequenceWithMods('unimod'), self.isoforms_seqmods, "Isoform4 not in the isoforms list")    
        self.assertNotIn(self.notisoform4.getSequenceWithMods('unimod'), self.isoforms_seqmods, "This should not be in the isoforms list")

    def test_calUIS(self) :
        self.theOtherIsoforms = [ self.isoform4_1 , self.isoform4_3 , self.isoform4_4 ]
        self.UISmass , self.UISannot = self.isoform4_2.cal_UIS(self.theOtherIsoforms, UISorder = 2,  ionseries = ['y'], fragmentlossgains = [0,], precision = 1e-5, frg_z_list = [1], mass_limits = [300,1500])
        self.assertEqual(len(self.UISmass), 7,  "wrong number of UIS. They should be 7!")
        self.assertAlmostIn((646.3406318939999,928.365933736), self.UISmass, "The UIS (646 , 928) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,1373.5984528780002), self.UISmass, "The UIS (646 , 1373) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,1373.5984528780002), self.UISmass, "The UIS (646 , 1373) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,813.338990706), self.UISmass, "The UIS (646 , 813) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,1316.5769891520001), self.UISmass, "The UIS (646 , 1316) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,1169.508575234), self.UISmass, "The UIS (646 , 1169) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,1098.471461444), self.UISmass, "The UIS (646 , 1098) is not in the UIS list!")
        self.assertAlmostIn((646.3406318939999,1041.449997718), self.UISmass, "The UIS (646 , 1041) is not in the UIS list!")
        
        
    def test_comparePeptideFragments(self):
        self.matched , self.unmatched = self.isoform1.comparePeptideFragments([self.isoform2], ['y','b'], precision = 1e-5)
        self.isoform1_b3_1 = self.isoform1.getMZfragment('b', 3, 1)
        self.isoform1_b7_2 = self.isoform1.getMZfragment('b', 7, 2)
        self.matchedMasses = [ mass for (_,_,_,_,mass) in self.matched ]
        self.unmatchedMasses = [ mass for (_,_,_,_,mass) in self.unmatched ]
        self.assertIn(self.isoform1_b3_1, self.matchedMasses, "ion b3_1 must be in the matched list!")
        self.assertNotIn(self.isoform1_b3_1, self.unmatchedMasses, "ion b3_1 must NOT be in the unmatched list!")
        self.assertNotIn(self.isoform1_b7_2, self.matchedMasses, "ion b7_2 must NOT be in the matched list!")
        self.assertIn(self.isoform1_b7_2, self.unmatchedMasses, "ion b3_1 must be in the unmatched list!")
        

    def test_getMassFromSequence(self):
        # TODO
        pass


    def test_getDeltaMassFromSequence(self):
        # TODO
        pass

    def test_pseudoreverse(self):
        # test that we dont loose any elements
        initial_len = len(self.mypep.sequence)
        initial_seq = self.mypep.sequence
        self.mypep.pseudoreverse()
        self.assertEqual( len(self.mypep.sequence), initial_len) 
        self.assertNotEqual( self.mypep.sequence, initial_seq) 

    def test_pseudoreverse_2(self):
        # test that we dont loose any elements
        initial_len = len(self.mypep.sequence)
        initial_seq = self.mypep.sequence
        newseq = self.mypep.pseudoreverse(initial_seq)
        self.assertEqual( len(newseq), initial_len) 
        self.assertNotEqual( newseq, initial_seq) 

    def test_shuffle_sequence(self):
        # test that we dont loose any elements
        initial_len = len(self.mypep.sequence)
        self.mypep.shuffle_sequence()
        self.assertEqual( len(self.mypep.sequence), initial_len) 

    def test_shuffle_sequence_1(self):
        # test that we also shuffle the modifications
        # TODO
        pass

    def test_get_decoy_Q3(self):
        # def test_get_decoy_Q3(self, frg_serie, frg_nr, frg_z, blackList=[], max_tries=1000):
        # TODO
        pass

    def test_getComposition(self):
        # TODO
        pass

    def test_getCompositionSeq(self):
        # def _getCompositionSeq(self, sequence, modifications = {}):
        # TODO
        pass

    def test_getAminoacidList(self):
        #def _getAminoacidList(self, fullList=False):
        pass
        print(self.mypep._getAminoacidList())

    def test_getMZ(self):
        self.assertAlmostEqual(self.mypep.getMZ(1), 1145.63470545)
        self.assertAlmostEqual(self.mypep.getMZ(2), 573.320990973)

        # TODO 'N15' 'SILAC_K6R10' 'SILAC_K8R10' 'SILAC_K8R6' :

        # TODO calculate true values
        # self.assertAlmostEqual(self.mypep.getMZ(1, "N15"), 1145.63470545)
        # self.assertAlmostEqual(self.mypep.getMZ(2, "N15"), 573.320990973)

    def test_getMZfragment(self):
        # def getMZfragment(self,ion_type,ion_number,ion_charge, label = '', fragmentlossgain = 0.0):
        self.assertAlmostEqual( self.mypep.getMZfragment("y", 4, 1), 478.2442291)
        self.assertAlmostEqual( self.mypep.getMZfragment("y", 5, 1), 577.312643018)
        self.assertAlmostEqual( self.mypep.getMZfragment("y", 6, 1), 664.344671428)
        self.assertRaises(SystemExit, self.mypep.getMZfragment, "y", 6, 1, "dummylabel")

        # TODO 'N15' 'SILAC_K6R10' 'SILAC_K8R10' 'SILAC_K8R6' :


    def test_test(self):
        peptide.test()

if __name__ == '__main__':
    unittest.main()
