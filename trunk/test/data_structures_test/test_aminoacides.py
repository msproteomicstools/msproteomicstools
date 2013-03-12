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

import msproteomicstoolslib.data_structures.aminoacides as aminoacides

class TestUnitAminoAcides(unittest.TestCase):

    def setUp(self):
        self.myAAs = aminoacides.Aminoacides()

    def test_getAminoacid(self):
        self.assertEqual(self.myAAs.getAminoacid("A").code3, "Ala")
        self.assertEqual(self.myAAs.getAminoacid("A").composition, {'H': 5, 'C': 3, 'O': 1, 'N': 1}  )
        self.assertAlmostEqual(self.myAAs.getAminoacid("A").deltaMass, 71.03711379)

    def test_getAminoacid_fail(self):
        self.assertRaises(Exception, self.myAAs.getAminoacid,"B")

    def test_getAminoacid_all(self):

        ##1-letter code 3-letter code Chemical formula Monoisotopic
        ##-----------------------------------------------------------------
        ##A           Ala              H5C3ON               71.03711379
        ##R           Arg              H12C6ON4             156.101111044
        ##N           Asn              H6C4O2N2             114.042927452
        ##D           Asp              H5C4O3N              115.02694303
        ##C           Cys              H5C3SON              103.00918448
        ##E           Glu              H7C5O3N              129.042593094
        ##Q           Gln              H8C5O2N2             128.058577516
        ##G           Gly              H3C2ON               57.021463726
        ##H           His              H7C6ON3              137.058911874
        ##I           Ile              H11C6ON              113.084063982
        ##L           Leu              H11C6ON              113.084063982
        ##K           Lys              H12C6ON2             128.094963024
        ##M           Met              H9C5SON              131.040484608
        ##F           Phe              H9C9ON               147.068413918
        ##P           Pro              H7C5ON               97.052763854
        ##S           Ser              H5C3O2N              87.03202841
        ##T           Thr              H7C4O2N              101.047678474
        ##W           Trp              H10C11ON2            186.07931296
        ##Y           tyr              H9C9O2N              163.063328538
        ##V           Val              H9C5ON               99.068413918  

        codes = "A R N D C E Q G H I L K M F P S T W Y V".split()
        codes3 = "Ala Arg Asn Asp Cys Glu Gln Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp tyr Val".split()
        chemical_comp = """H5C3ON   
            H12C6ON4 
            H6C4O2N2 
            H5C4O3N  
            H5C3SON
            H7C5O3N  
            H8C5O2N2 
            H3C2ON   
            H7C6ON3  
            H11C6ON  
            H11C6ON  
            H12C6ON2 
            H9C5SON  
            H9C9ON   
            H7C5ON   
            H5C3O2N  
            H7C4O2N  
            H10C11ON2
            H9C9O2N  
            H9C5ON   """.split()

        composition =  [
        
            {'H': 5, 'C': 3, 'O': 1, 'N': 1},
            {'H': 12, 'C': 6, 'O': 1, 'N': 4},
            {'H': 6, 'C': 4, 'O': 2, 'N': 2},
            {'H': 5, 'C': 4, 'O': 3, 'N': 1},
            {'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1},
            {'H': 7, 'C': 5, 'O': 3, 'N': 1},
            {'H': 8, 'C': 5, 'O': 2, 'N': 2},
            {'H': 3, 'C': 2, 'O': 1, 'N': 1},
            {'H': 7, 'C': 6, 'O': 1, 'N': 3},
            {'H': 11, 'C': 6, 'O': 1, 'N': 1},
            {'H': 11, 'C': 6, 'O': 1, 'N': 1},
            {'H': 12, 'C': 6, 'O': 1, 'N': 2},
            {'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1},
            {'H': 9, 'C': 9, 'O': 1, 'N': 1},
            {'H': 7, 'C': 5, 'O': 1, 'N': 1},
            {'H': 5, 'C': 3, 'O': 2, 'N': 1},
            {'H': 7, 'C': 4, 'O': 2, 'N': 1},
            {'H': 10, 'C': 11, 'O': 1, 'N': 2},
            {'H': 9, 'C': 9, 'O': 2, 'N': 1},
            {'H': 9, 'C': 5, 'O': 1, 'N': 1}

        ]

        monoisotopic_mass = [
            71.03711379,
            156.101111044 ,
            114.042927452 ,
            115.02694303  ,
            103.00918448  ,
            129.042593094 ,
            128.058577516 ,
            57.021463726  ,
            137.058911874 ,
            113.084063982 ,
            113.084063982 ,
            128.094963024 ,
            131.040484608 ,
            147.068413918 ,
            97.052763854  ,
            87.03202841   ,
            101.047678474 ,
            186.07931296  ,
            163.063328538 ,
            99.068413918 ]

        for i,code in enumerate(codes):
            aa = self.myAAs.getAminoacid(code)
            # print "%-10s  %-15s  %-20s %s %s" % (aa.code, aa.code3, aa.compositionString(), aa.deltaMass, i)
            self.assertEqual(aa.code3, codes3[i])
            self.assertEqual(aa.composition, composition[i])
            self.assertAlmostEqual(aa.deltaMass, monoisotopic_mass[i])


if __name__ == '__main__':
    unittest.main()
