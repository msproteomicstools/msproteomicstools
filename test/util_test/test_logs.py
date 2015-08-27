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
import unittest
import os

from msproteomicstoolslib.format.ProteinDB import ProteinDB
import msproteomicstoolslib.util.logs


class multiprocessStuff():
    def __init__(self, logger):
        self.logger = logger


class TestMultiProcessingLog(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.fastafile = os.path.join(self.dirname, "..", "data", 'smallDB.fasta')
        self.db = ProteinDB()
        self.db.readFasta(self.fastafile)
        self.P31946 = self.db.proteinDictionary['P31946']

    def test_readFasta(self):
        self.code1 = self.P31946.code1
        self.code2 = self.P31946.code2
        self.modres = self.P31946.modres
        self.ncbi_tax_id = self.P31946.ncbi_tax_id
        self.description = self.P31946.description
        self.sequence = self.P31946.sequence

        print(self.code1)
        print(self.code2)
        print(self.modres)
        print(self.ncbi_tax_id)
        print(self.description)
        print(self.sequence)


    def test_get_proteins_containing_peptide(
            self):  #To-Do : test with some peptides present and not present. Also shared peptides.
        pass

    def test_pseudoreverseDB(
            self):  #To-Do : Pseudo-reverse and check a couple of reversed proteins. Pseudo-reverse the DB twice, and check the DB is equal.
        pass


class TestUnitProtein(unittest.TestCase):
    def setUp(self):
        #>sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens GN=YWHAB PE=1 SV=3
        #MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS
        #WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY
        #LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY
        #YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD
        #AGEGEN
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.fastafile = os.path.join(self.dirname, "..", "data", 'smallDB.fasta')
        self.db = ProteinDB()
        self.db.readFasta(self.fastafile)
        self.P31946 = self.db.proteinDictionary['P31946']
        self.trypsin = {'terminus': 'C', 'cleave': ['K', 'R'], 'exceptions': ['KP', 'RP']}
        self.Lys_N = {'terminus': 'N', 'cleave': ['K'], 'exceptions': []}
        self.pep_tryp1 = 'TAFDEAIAELDTLNEESYK'
        self.pep_Lys_N = 'KGDYFRYLSEVASGDN'


    def test_proteinWeight(self):
        pass


    def test_digest(self):
        self.P31946_peptides_trypsin = self.P31946.digest(self.trypsin)

        self.assertIn(self.pep_tryp1, self.P31946_peptides_trypsin,
                      "A tryptic peptide has not been found after tryptic digestion of a protein!")
        self.assertNotIn(self.pep_Lys_N, self.P31946_peptides_trypsin,
                         "A non tryptic peptide has been found after tryptic digestion of a protein!")

        self.P31946_peptides_Lys_N = self.P31946.digest(self.Lys_N)
        self.assertNotIn(self.pep_tryp1, self.P31946_peptides_Lys_N,
                         "A  tryptic peptide has been found after tryptic digestion of a protein!")
        self.assertIn(self.pep_Lys_N, self.P31946_peptides_Lys_N,
                      "A Lys-N peptide has not been found after tryptic digestion of a protein!")


if __name__ == '__main__':
    unittest.main()
