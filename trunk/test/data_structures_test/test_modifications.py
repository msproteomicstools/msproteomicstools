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

import msproteomicstoolslib.data_structures.modifications as modifications

def fail():
        raise ValueError('Misspellled errrorr messageee')

class TestUnitModifications(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        joined = os.path.join(dirname, "..")
        joined = os.path.join(joined, "data")
        self.filename = os.path.join(joined, "modsfile_test.txt")
        self.mods = modifications.Modifications()

    def test_readModificationsFile(self) :
        self.mods.readModificationsFile(self.filename)

    def test_translateModificationsFromSequence(self):
        self.mods.translateModificationsFromSequence("PEPTIDEMEK", "TPP")
        self.mods.translateModificationsFromSequence("PEPTIDEM[147]EK", "TPP")
        self.mods.translateModificationsFromSequence("PEPTIDEM(UniMod:35)EK", "unimod")

        # should nor work with dummy input
        self.assertRaises(SystemExit, self.mods.translateModificationsFromSequence,"P", "Dummy")
        self.assertRaises(SystemExit, self.mods.translateModificationsFromSequence,"PEPTIDEM(unimod:11111)K", "unimod")
        self.assertRaises(SystemExit, self.mods.translateModificationsFromSequence,"n[11111]PEPTIDEMK", "TPP")

        self.mods.translateModificationsFromSequence("n[43]PEPTIDEM[147]K", "TPP")
        self.mods.translateModificationsFromSequence("(UniMod:5)PEPTIDEM(UniMod:35)K", "unimod")

    def test_printModifications(self):
        self.mods.printModifications()

    def test_Modificationtest(self):
        modifications.test([])

class TestUnitModification_2(unittest.TestCase):

    def setUp(self):
        self.mod = modifications.Modification('R', 'R[166]', 267, '[+10]', True, {'C' : -6 , '13C' : 6 , 'N' : -4 , '15N' : 4 } )  

    def test_getcode(self):
        self.assertEqual( self.mod.getcode("TPP"), 'R[166]')
        self.assertEqual( self.mod.getcode("unimod"), 'R(UniMod:267)')
        self.assertEqual( self.mod.getcode("ProteinPilot"), 'R[+10]')
        self.assertRaises(SystemExit, self.mod.getcode,"Dummy")

if __name__ == '__main__':
    unittest.main()
