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

import msproteomicstoolslib.data_structures.elements as elements

class TestUnitElement(unittest.TestCase):

    def setUp(self):
        self.isots = elements.Elements()

    def old_function_printout(self):
        # TODO what do we test here?
        isots = elements.Elements()
        
        for el in isots.list:
            print(el.symbol)
            totalAb = 0
            for k,m in zip(el.isotMass,el.isotAbundance):
                print("Mass: %s, Abundance: %s" % (k,m))
                totalAb += m
            print("Total abundance : " , totalAb)
                
        #Monoisotopic masses
        for el in isots.list:
            print(el.symbol , el.isotMass[0])

    def test_getElement(self):
        self.assertAlmostEqual(self.isots.getElement("C").isotMass[0], 12)
        self.assertAlmostEqual(self.isots.getElement("H").isotMass[0], 1.007825032)
        self.assertAlmostEqual(self.isots.getElement("N").isotMass[0], 14.00307401)
        self.assertAlmostEqual(self.isots.getElement("O").isotMass[0], 15.99491462)
        self.assertAlmostEqual(self.isots.getElement("S").isotMass[0], 31.97207069)
        self.assertAlmostEqual(self.isots.getElement("P").isotMass[0], 30.97376151)

        self.assertRaises(Exception, self.isots.getElement,"x") 

    def test_addElement(self):
      # def addElement(self,symbol,isotMass,isotAbundance):
      self.isots.addElement("C", [12], [1.0])

    def test_test(self):
        elements.test()

class TestUnitFormulas(unittest.TestCase):

    def setUp(self):
        self.formulas = elements.Formulas()

    def test_mass(self):
        self.assertAlmostEqual(self.formulas.mass({ "H" : 2}),  2.015650064)
        self.assertAlmostEqual(self.formulas.mass({ "C13" : 1}),  13.00335484)
        self.assertRaises(SystemExit, self.formulas.mass,{ "C18" : 1})
        self.assertRaises(SystemExit, self.formulas.mass,{ "C18C45" : 1})

    def test_add2components(self):

        from msproteomicstoolslib.data_structures.elements import Formulas
        self.assertEqual( {'H': 5, 'O': 5, 'P': 1}, Formulas.add2components(Formulas.H2O, Formulas.H3PO4))
        self.assertEqual( {'H': 1, 'O': 3, 'P': 1}, Formulas.substract2components(Formulas.H3PO4, Formulas.H2O) )
        self.assertEqual( 114.01844984 , Formulas.mass({'C13' : 6 , '18O' : 2}) )

    def test_compositionString(self):
        self.assertEqual( 'H5O5P', elements.Formulas.compositionString( {'H': 5, 'O': 5, 'P': 1}) )



if __name__ == '__main__':
    unittest.main()
