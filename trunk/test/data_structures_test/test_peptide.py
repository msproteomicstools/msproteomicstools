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

import unittest
import os

import msproteomicstoolslib.data_structures.peptide as peptide
# from msproteomicstoolslib.data_structures.peptide import *

class TestUnitPeptide(unittest.TestCase):

    def setUp(self):
	    self.mypep = peptide.Peptide('LIGPTSVVMGR',modifications = {9:15.99491462})
	    self.mypep2 = peptide.Peptide('LIGPTSVVMGR')

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
		## 	masses.append(str(v))
		## 	for ionserie in pep.iontypes :
		## 		masses.append(str(pep.getMZfragment(ionserie,v,1)))
		## 	print ", ".join(masses)
		## 	masses = []
		## 	
		## 	#for property, value in vars(pep).iteritems() :
		## 	#	print property , " : " , value

    def test_peptide_with_modifications(self):

		# check that the total mass, mz and composition is correct
		pep = peptide.Peptide('LIGPTSVVMGR',modifications = {9:15.99491462})
		self.assertAlmostEqual(pep.mass, 1144.62742895)
		self.assertEqual(pep._getComposition(), {'H': 86, 'C': 49, 'S': 1, 'O': 13, 'N': 14} )
		self.assertAlmostEqual(pep.getMZ(1), 1145.63470545)
		
		# TODO also check some ion series
		## print '#,' , ", ".join(pep.iontypes)
		## masses = []
		## for v in range(1,len(pep.sequence)+1) :
		## 	masses.append(str(v))
		## 	for ionserie in pep.iontypes :
		## 		masses.append(str(pep.getMZfragment(ionserie,v,1)))
		## 	print ", ".join(masses)
		## 	masses = []
		## 	
		## 	#for property, value in vars(pep).iteritems() :
		## 	#	print property , " : " , value

    def test_getMZfragment(self):
        # def getMZfragment(self,ion_type,ion_number,ion_charge, label = '', fragmentlossgain = 0.0):
        self.assertAlmostEqual( self.mypep.getMZfragment("y", 4, 1), 478.2442291)
        self.assertAlmostEqual( self.mypep.getMZfragment("y", 5, 1), 577.312643018)
        self.assertAlmostEqual( self.mypep.getMZfragment("y", 6, 1), 664.344671428)

    def test_getMZ(self):
		self.assertAlmostEqual(self.mypep.getMZ(1), 1145.63470545)
		self.assertAlmostEqual(self.mypep.getMZ(2), 573.320990973)

        # TODO calculate true values
		# self.assertAlmostEqual(self.mypep.getMZ(1, "N15"), 1145.63470545)
		# self.assertAlmostEqual(self.mypep.getMZ(2, "N15"), 573.320990973)

    def test_shuffle_sequence(self):
        # test that we dont loose any elements
        initial_len = len(self.mypep.sequence)
        self.mypep.shuffle_sequence()
        self.assertEqual( len(self.mypep.sequence), initial_len) 

    def test_shuffle_sequence_1(self):
        # test that we also shuffle the modifications
        pass

    def test_pseudoreverse(self):
        # test that we dont loose any elements
        initial_len = len(self.mypep.sequence)
        initial_seq = self.mypep.sequence
        self.mypep.pseudoreverse()
        self.assertEqual( len(self.mypep.sequence), initial_len) 
        self.assertNotEqual( self.mypep.sequence, initial_seq) 

    def test_getDeltaMassFromSequence(self):
        # TODO
        pass

if __name__ == '__main__':
    unittest.main()
