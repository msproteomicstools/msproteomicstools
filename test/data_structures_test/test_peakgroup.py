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

import msproteomicstoolslib.data_structures.PeakGroup as peakgroup

class MockRun():

    def __init__(self):
        pass

class MockPeptide():

    def __init__(self):
        pass

    def setClusterID(self, feature_id, id_):
        pass

class TestUnitPeakGroupBase(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_pg(self):
        pg = peakgroup.PeakGroupBase()
        self.assertTrue(True)

    def test_set_fdr_score(self):
        pg = peakgroup.PeakGroupBase()
        pg.set_fdr_score(0.1)
        self.assertAlmostEqual(pg.get_fdr_score(), 0.1)

    def test_set_normalized_retentiontime(self):
        pg = peakgroup.PeakGroupBase()
        pg.set_normalized_retentiontime(0.1)
        self.assertAlmostEqual(pg.get_normalized_retentiontime(), 0.1)

    def test_set_feature_id(self):
        pg = peakgroup.PeakGroupBase()
        pg.set_feature_id("f_1")
        self.assertEqual(pg.get_feature_id(), "f_1")

    def test_set_intensity(self):
        pg = peakgroup.PeakGroupBase()
        pg.set_intensity(0.1)
        self.assertAlmostEqual(pg.get_intensity(), 0.1)

    def test_select(self):
        pg = peakgroup.PeakGroupBase()
        pg.select_this_peakgroup()
        self.assertTrue(True)

class TestUnitMinimalPeakGroup(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_pg(self):
        pg = peakgroup.MinimalPeakGroup("f_1", 0.1, 100, False, -1, None)
        self.assertTrue(True)

    def test_setClusterId(self):
        peptide = MockPeptide()
        pg = peakgroup.MinimalPeakGroup("f_1", 0.1, 100, False, -1, peptide)
        pg.setClusterID(5)
        self.assertEqual(pg.get_cluster_id(), 5)

    def test_exception(self):
        peptide = MockPeptide()
        pg = peakgroup.MinimalPeakGroup("f_1", 0.1, 100, False, -1, peptide)

        self.assertRaises(Exception, pg.set_fdr_score, 4)
        self.assertRaises(Exception, pg.set_normalized_retentiontime, 4)
        self.assertRaises(Exception, pg.set_feature_id, 4)
        self.assertRaises(Exception, pg.set_intensity, 4)

    def test_dscore(self):
        peptide = MockPeptide()
        pg = peakgroup.MinimalPeakGroup("f_1", 0.1, 100, False, -1, peptide, None, 42)
        self.assertAlmostEqual(pg.get_dscore(), 42)

class TestUnitGuiPeakGroup(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_pg(self):
        pg = peakgroup.GuiPeakGroup(0.1, 100, 50, 60, 7, None)

    def test_get_value(self):
        pg = peakgroup.GuiPeakGroup(0.1, 100, 50, 60, 8.0, None)
        self.assertAlmostEqual(pg.get_value("m_score"), 0.1)
        self.assertAlmostEqual(pg.get_value("Intensity"), 100)
        self.assertAlmostEqual(pg.get_value("leftWidth"), 50)
        self.assertAlmostEqual(pg.get_value("rightWidth"), 60)
        self.assertRaises(Exception, pg.get_value, "dummy")

class TestUnitGeneralPeakGroup(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_pg(self):
        pg = peakgroup.GeneralPeakGroup([], None, None)

    def test_get_value(self):
        run = MockRun()
        run.header_dict = {"dummy": 0, "d_score": 1}
        pg = peakgroup.GeneralPeakGroup([ 8, 9.0 ], run, None)
        pg.set_value("dummy", 42)
        self.assertAlmostEqual(pg.get_value("dummy"), 42)
        pg.set_value("dummy", None)
        self.assertEqual(pg.get_value("dummy"), "NA")
        self.assertAlmostEqual(pg.get_dscore(), 9.0)

    def test_setClusterId(self):
        pg = peakgroup.GeneralPeakGroup([], None, None)
        pg.setClusterID(5)
        self.assertEqual(pg.get_cluster_id(), 5)

if __name__ == '__main__':
    unittest.main()
