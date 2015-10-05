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
import os

import unittest
from nose.plugins.attrib import attr

import msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm as algo
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
import msproteomicstoolslib.data_structures.Precursor as precursor
import msproteomicstoolslib.data_structures.PrecursorGroup as precursor_group

class MockRun():

    def __init__(self, id_):
        self.id_ = id_
        self.orig_filename = "test"

    def get_id(self):
        return self.id_


class TestUnitAlignmentAlgo(unittest.TestCase):

    def setUp(self):

        # 0. id
        # 1. quality score (FDR)
        # 2. retention time (normalized)
        # 3. intensity

        m = Multipeptide()
        m.set_nr_runs(2)

        # Run 1
        r = MockRun("0_1")
        p = precursor.Precursor("precursor_1", r)
        pg_tuple = ("someID_1", 0.1, 100, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), r)
        prgr.addPrecursor(p)
        m.insert("0_1", prgr)

        # Run 2:
        #  - peakgroup 2 : RT = 105 seconds
        #  - peakgroup 3 : RT = 120 seconds
        r = MockRun("0_2")
        p = precursor.Precursor("precursor_1", r)
        pg_tuple = ("peakgroup2", 0.2, 105, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        pg_tuple = ("peakgroup3", 0.18, 130, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), r)
        prgr.addPrecursor(p)
        m.insert("0_2", prgr)

        self.mpep = m
        self.al = algo.AlignmentAlgorithm()
        self.al.verbose = True

    def test_cluster_1(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 10
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.15

        self.al._align_features_cluster(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 1)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )

    def test_cluster_twoPG(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 10
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.3

        self.al._align_features_cluster(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 2)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertEqual( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup().get_feature_id(), "peakgroup2")

    def test_cluster_twoPG_large(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 40
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.3

        self.al._align_features_cluster(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 2)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertEqual( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup().get_feature_id(), "someID_1")
        self.assertEqual( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup().get_feature_id(), "peakgroup3")

    def test_cluster_toplvl(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 10
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.15

        self.al.align_features([ self.mpep ], rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "best_cluster_score")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 1)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )

    def test_best_1(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 10
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.15

        self.al._align_features_best(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 1)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )

    def test_best_twoPG(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 10
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.3

        self.al._align_features_best(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 2)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertEqual( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup().get_feature_id(), "peakgroup2")

    def test_best_twoPG_large(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 40
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.3

        self.al._align_features_best(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 2)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertEqual( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup().get_feature_id(), "peakgroup3")

    def test_best_3(self):
        """Test the best overall algorithm"""

        aligned_fdr_cutoff = 0.3
        fdr_cutoff = 0.15
        rt_diff_cutoff = 4

        self.al._align_features_best(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "dummy")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 1)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )

    def test_best_4(self):
        """Test the best overall algorithm"""

        aligned_fdr_cutoff = 0.3
        fdr_cutoff = 0.15
        rt_diff_cutoff = 4

        self.al._align_features_best(self.mpep, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "best_overall")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 1)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )

    def test_best_toplvl(self):
        """Test the best overall algorithm"""

        rt_diff_cutoff = 10
        fdr_cutoff = 0.15
        aligned_fdr_cutoff = 0.15

        self.al.align_features([ self.mpep ], rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, "best_overall")
        self.assertEqual( len(self.mpep.get_selected_peakgroups()), 1)
        self.assertIsNotNone( self.mpep.getPrecursorGroup("0_1").getPrecursor("precursor_1").get_selected_peakgroup() )
        self.assertIsNone( self.mpep.getPrecursorGroup("0_2").getPrecursor("precursor_1").get_selected_peakgroup() )

if __name__ == '__main__':
    unittest.main()
