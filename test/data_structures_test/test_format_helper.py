#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013-2017, ETH Zurich
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

import msproteomicstoolslib.data_structures.FormatHelper as FormatHelper

class TestUnitFormatHelper(unittest.TestCase):

    def setUp(self):
        self.h = FormatHelper()

        self.test_k1 = "DECOY_44736_NVEVIEDDKQGIIR/2_y12"
        self.test_k2 = "DECOY_155153_GYEDPPAALFR/2_y7_2"
        self.test_k3 = "DECOY_44736_y6_1_NVEVIEDDKQGIIR_2"
        self.test_k4 = "1002781_TGLC(UniMod:4)QFEDAFTQLSGATPIGAGIDAR_3"

    def test_check_format(self):
        self.assertTrue(self.h._has_openswath_format(self.test_k1))
        self.assertTrue(self.h._has_openswath_format(self.test_k2))
        self.assertTrue(self.h._has_openswath_format(self.test_k3))
        self.assertTrue(self.h._has_openswath_format(self.test_k4))

    def test_format_parse(self):
        # returns tuple (decoy, trgr_nr, sequence, prec_charge, fr_id, fr_charge)

        p1 = self.h.parse( self.test_k1 )
        p2 = self.h.parse( self.test_k2 )
        p3 = self.h.parse( self.test_k3 )
        p4 = self.h.parse( self.test_k4 )

        self.assertTrue(p1 is not None)
        self.assertTrue(p2 is not None)
        self.assertTrue(p3 is not None)
        self.assertTrue(p4 is not None)

        self.assertTrue(p1[0])
        self.assertTrue(p2[0])
        self.assertTrue(p3[0])
        self.assertTrue(not p4[0])


        self.assertTrue(p1[1] == "44736")
        self.assertTrue(p2[1] == "155153")
        self.assertTrue(p3[1] == "44736")
        self.assertTrue(p4[1] == "1002781")

        self.assertTrue(p1[2] == "NVEVIEDDKQGIIR")
        self.assertTrue(p2[2] == "GYEDPPAALFR")
        self.assertTrue(p3[2] == "NVEVIEDDKQGIIR")
        self.assertTrue(p4[2] == "TGLC(UniMod:4)QFEDAFTQLSGATPIGAGIDAR")

        self.assertTrue(p1[3] == "2")
        self.assertTrue(p2[3] == "2")
        self.assertTrue(p3[3] == "2")
        self.assertTrue(p4[3] == "3")

    def test_compute_transitiongroup_from_key(self):
        print (self.h._compute_transitiongroup_from_key(self.test_k1))
        self.assertEqual(self.h._compute_transitiongroup_from_key(self.test_k1), "DECOY_NVEVIEDDKQGIIR/2")
        self.assertEqual(self.h._compute_transitiongroup_from_key(self.test_k2), "DECOY_GYEDPPAALFR/2")
        self.assertEqual(self.h._compute_transitiongroup_from_key(self.test_k3), "DECOY_NVEVIEDDKQGIIR/2")
        self.assertEqual(self.h._compute_transitiongroup_from_key(self.test_k4), "TGLC(UniMod:4)QFEDAFTQLSGATPIGAGIDAR/3")

if __name__ == '__main__':
    unittest.main()
