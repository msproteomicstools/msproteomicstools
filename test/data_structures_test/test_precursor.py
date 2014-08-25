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

import msproteomicstoolslib.data_structures.Precursor as precursor
from msproteomicstoolslib.data_structures.PeakGroup import PeakGroupBase

class MockPeakGroup():

    def __init__(self, fdr=1, id=1, selected=False):
        self.fdr_score = fdr
        self.id = id
        self.sel = selected

    def get_id(self):
        return self.id

    def is_selected(self):
        return self.id


class TestUnitPrecursorBase(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_precursor(self):
        self.assertRaises(Exception, precursor.PrecursorBase, 0, 0)

class TestUnitPrecursor(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_precursor(self):
        p = precursor.Precursor("precursor_2", [])
        str(p)
        self.assertTrue(True)

    def test_get_id(self):
        p = precursor.Precursor("precursor_2", [])
        self.assertEqual(p.get_id(), "precursor_2")

    def test_get_decoy(self):
        p = precursor.Precursor("precursor_2", [])
        self.assertFalse(p.get_decoy())
        p.set_decoy("TRUE")
        self.assertTrue(p.get_decoy())
        p.set_decoy("FALSE")
        self.assertFalse(p.get_decoy())
        p.set_decoy("1")
        self.assertTrue(p.get_decoy())
        p.set_decoy("0")
        self.assertFalse(p.get_decoy())

        self.assertRaises(Exception, p.set_decoy, "dummy")

    def test_add_peakgroup_tpl(self):
        """
        0. id
        1. quality score (FDR)
        2. retention time (normalized)
        3. intensity
        (4. d_score optional)
        """

        pg_tuple = ("someID", 0.1, 100, 10000, 2)

        p = precursor.Precursor("precursor_2", [])
        self.assertRaises(Exception, p.add_peakgroup_tpl, pg_tuple, "notMatchingID", 4)
        p.add_peakgroup_tpl(pg_tuple, "precursor_2", 4)

        self.assertEqual( len(list(p.get_all_peakgroups())), 1)
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), 4)
        self.assertAlmostEqual( firstpg.get_fdr_score(), 0.1)
        self.assertAlmostEqual( firstpg.get_normalized_retentiontime(), 100)
        self.assertAlmostEqual( firstpg.get_intensity(), 10000)

        pg_tuple = ("someID", 0.1, 100, 10000)
        p = precursor.Precursor("precursor_2", [])
        p.add_peakgroup_tpl(pg_tuple, "precursor_2", 4)

        self.assertEqual( len(list(p.get_all_peakgroups())), 1)

    def test_setClusterId(self):
        pg_tuple = ("someID", 0.1, 100, 10000, 2)
        p = precursor.Precursor("precursor_2", [])
        p.add_peakgroup_tpl(pg_tuple, "precursor_2", 4)
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), 4)

        p.setClusterID("someID", 5)
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), 5)

    def test_selectpg(self):
        pg_tuple = ("someID", 0.1, 100, 10000, 2)
        p = precursor.Precursor("precursor_2", [])
        p.add_peakgroup_tpl(pg_tuple, "precursor_2", 4)
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), 4)

        p.unselect_pg("someID")
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), -1)
        p.select_pg("someID")
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), 1)
        p.unselect_pg("someID")
        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), -1)

    def test_selection(self):
        p = precursor.Precursor("precursor_2", [])
        self.assertIsNone( p.get_best_peakgroup() )

        pg_tuple = ("someID", 0.1, 100, 10000, 2)
        p.add_peakgroup_tpl(pg_tuple, "precursor_2", 1)
        pg_tuple = ("someID_", 0.01, 105, 10000, 2)
        p.add_peakgroup_tpl(pg_tuple, "precursor_2", -1)

        self.assertEqual( p.get_selected_peakgroup().get_feature_id(), "someID")
        self.assertEqual( p.get_best_peakgroup().get_feature_id(), "someID_")
        self.assertEqual( p.find_closest_in_iRT(99).get_feature_id(), "someID")
        self.assertEqual( p.find_closest_in_iRT(110).get_feature_id(), "someID_")
        self.assertEqual( p.find_closest_in_iRT(102.8).get_feature_id(), "someID_")
        self.assertEqual( len(list(p.getClusteredPeakgroups())), 1)
        self.assertEqual( list(p.getClusteredPeakgroups())[0].get_feature_id(), "someID")

        # Un-select all pg
        # Only one pg should be selected at a time
        p.unselect_all()
        self.assertIsNone(p.get_selected_peakgroup())
        p.select_pg("someID_")
        self.assertRaises(AssertionError, p.select_pg, "someID")

class TestUnitGeneralPrecursor(unittest.TestCase):

    def setUp(self):
        pass
        
    def test_create_precursor(self):
        p = precursor.GeneralPrecursor("precursor_2", [])
        self.assertTrue(True)

    def test_get_id(self):
        p = precursor.GeneralPrecursor("precursor_2", [])
        self.assertEqual(p.get_id(), "precursor_2")

    def test_append(self):

        # TODO is this even used anywhere ??? 
        pg = MockPeakGroup(1, "precursor_2")
        p = precursor.GeneralPrecursor("precursor_2", [])
        p.append(pg)
        self.assertEqual( len(list(p.get_all_peakgroups())), 1)

    def test_add_peakgroup(self):

        pg = PeakGroupBase()
        p = precursor.GeneralPrecursor("precursor_2", [])
        p.add_peakgroup(pg)
        self.assertEqual( len(list(p.get_all_peakgroups())), 1)

    def test_selectpg(self):
        pg = PeakGroupBase()
        pg.cluster_id_ = 4
        pg.id_ = "someID"
        p = precursor.GeneralPrecursor("precursor_2", [])
        p.add_peakgroup(pg)
        self.assertEqual( len(list(p.get_all_peakgroups())), 1)

        firstpg = list(p.get_all_peakgroups())[0]
        self.assertEqual( firstpg.get_cluster_id(), 4)

    def test_selection(self):
        p = precursor.GeneralPrecursor("precursor_2", [])
        self.assertIsNone( p.get_selected_peakgroup() )
        self.assertIsNone( p.get_best_peakgroup() )

        pg = PeakGroupBase()
        pg.cluster_id_ = 1
        pg.id_ = "someID"
        pg.normalized_retentiontime = 100
        p.add_peakgroup(pg)
        self.assertEqual( len(list(p.get_all_peakgroups())), 1)

        pg = PeakGroupBase()
        pg.cluster_id_ = 2
        pg.id_ = "someID_"
        pg.normalized_retentiontime = 105
        p.add_peakgroup(pg)
        self.assertEqual( len(list(p.get_all_peakgroups())), 2)

        self.assertEqual( p.get_selected_peakgroup().get_feature_id(), "someID")
        self.assertEqual( p.get_best_peakgroup().get_feature_id(), "someID_")
        self.assertEqual( p.find_closest_in_iRT(99).get_feature_id(), "someID")
        self.assertEqual( p.find_closest_in_iRT(110).get_feature_id(), "someID_")
        self.assertEqual( p.find_closest_in_iRT(102.8).get_feature_id(), "someID_")


if __name__ == '__main__':
    unittest.main()
