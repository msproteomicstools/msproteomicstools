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

import msproteomicstoolslib.algorithms.alignment.AlignmentMST as algo
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner

class MockRun():

    def __init__(self, id_):
        self.id_ = id_
        self.orig_filename = "test"

    def get_id(self):
        return self.id_

class Dummy():
    pass

class TestUnitAlignmentAlgo(unittest.TestCase):

    def setUp(self):

        import msproteomicstoolslib.data_structures.Precursor as precursor
        import msproteomicstoolslib.data_structures.PrecursorGroup as precursor_group
        import msproteomicstoolslib.format.TransformationCollection as transformations
        from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
        import msproteomicstoolslib.algorithms.alignment.AlignmentHelper as helper

        # 0. id
        # 1. quality score (FDR)
        # 2. retention time (normalized)
        # 3. intensity

        mpeps = [Multipeptide() for i in range(3)]
        [m.set_nr_runs(5) for m in mpeps]

        # Parameters
        self.initial_alignment_cutoff = 0.001

        runs = [MockRun("0_%s" % (i+1)) for i in range(5)]
        ids = 0
        for i in range(5):

            # Two alignment peptides
            p = precursor.Precursor("anchorpeptide_1", runs[i] )
            pg_tuple = ("id_%s" % ids, 0.0001, 100 + i*10, 10000)
            p.add_peakgroup_tpl(pg_tuple, "anchorpeptide_1", -1)
            prgr = precursor_group.PrecursorGroup(p.get_id(), runs[i])
            prgr.addPrecursor(p)
            mpeps[0].insert(runs[i].get_id(), prgr)
            ids += 1

            p = precursor.Precursor("anchorpeptide_2", runs[i] )
            pg_tuple = ("id_%s" % ids, 0.0001, 1000 + i*100, 10000)
            p.add_peakgroup_tpl(pg_tuple, "anchorpeptide_2", -1)
            prgr = precursor_group.PrecursorGroup(p.get_id(), runs[i])
            prgr.addPrecursor(p)
            mpeps[1].insert(runs[i].get_id(), prgr)
            ids += 1

            # The noise peptide
            p = precursor.Precursor("anchorpeptide_3", runs[i] )
            pg_tuple = ("id_%s" % ids, 0.0001, 500 + i*40, 10000)
            p.add_peakgroup_tpl(pg_tuple, "anchorpeptide_3", -1)
            prgr = precursor_group.PrecursorGroup(p.get_id(), runs[i])
            prgr.addPrecursor(p)
            mpeps[2].insert(runs[i].get_id(), prgr)
            ids += 1

        m = Multipeptide()
        m.set_nr_runs(5)

        # Run 1
        #  - peakgroup 1 : RT = 110 seconds [correct]
        p = precursor.Precursor("precursor_1", runs[0])
        pg_tuple = ("peakgroup1", 0.01, 100, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), runs[0])
        prgr.addPrecursor(p)
        m.insert(runs[0].get_id(), prgr)

        # Run 2:
        #  - peakgroup 2 : RT = 115 seconds [correct]
        #  - peakgroup 3 : RT = 130 seconds
        p = precursor.Precursor("precursor_1", runs[1])
        pg_tuple = ("peakgroup2", 0.2, 112, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        pg_tuple = ("peakgroup3", 0.18, 130, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), runs[1])
        prgr.addPrecursor(p)
        m.insert(runs[1].get_id(), prgr)

        # Run 3:
        #  - peakgroup 4 : RT = 120 seconds [correct]
        #  - peakgroup 5 : RT = 130 seconds
        p = precursor.Precursor("precursor_1", runs[2])
        pg_tuple = ("peakgroup4", 0.2, 120, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        pg_tuple = ("peakgroup5", 0.17, 130, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), runs[2])
        prgr.addPrecursor(p)
        m.insert(runs[2].get_id(), prgr)

        # Run 4:
        #  - peakgroup 6 : missing          [correct]
        #  - peakgroup 7 : RT = 145 seconds
        p = precursor.Precursor("precursor_1", runs[3])
        pg_tuple = ("peakgroup7", 0.18, 145, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), runs[3])
        prgr.addPrecursor(p)
        m.insert(runs[3].get_id(), prgr)

        # Run 5:
        #  - peakgroup 8 : RT = 140 seconds [correct]
        #  - peakgroup 9 : missing
        p = precursor.Precursor("precursor_1", runs[4])
        pg_tuple = ("peakgroup8", 0.1, 139, 10000)
        p.add_peakgroup_tpl(pg_tuple, "precursor_1", -1)
        prgr = precursor_group.PrecursorGroup(p.get_id(), runs[4])
        prgr.addPrecursor(p)
        m.insert(runs[4].get_id(), prgr)

        self.mpep = m
        self.exp = Dummy()
        self.exp.runs = runs

        mpeps.append(m)
        self.multipeptides = mpeps

        # Align all against all
        self.tr_data = transformations.LightTransformationData()
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        for run_0 in self.exp.runs:
            for run_1 in self.exp.runs:
                helper.addDataToTrafo(self.tr_data, run_0, run_1, spl_aligner, self.multipeptides, "linear", 30)

    def test_prepare(self):
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(algo.getDistanceMatrix(self.exp, self.multipeptides, spl_aligner))
        self.assertEqual(tree, [(3, 4), (2, 3), (1, 2), (0, 1)] )

    def test_alignBestCluster_0(self):
        """Test the best cluster align
        
        This is using the best possible conditions with only 7 seconds retention time cutoff

          - Run1 : 100s     [threadRT = 100s] 
          - Run2 : 112s     [threadRT = 106s]
          - Run3 : 120s     [threadRT = 118s]
          - Run4 : xxx      [threadRT = 126s]  (should be around 130s)
          - Run5 : 139s     [threadRT = 133s]
        """

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(algo.getDistanceMatrix(self.exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.exp.runs[a].get_id(), self.exp.runs[b].get_id()) for a,b in tree]

        alignment = algo.TreeConsensusAlignment(max_rt_diff = 6, fdr_cutoff = 0.1, aligned_fdr_cutoff = 0.25, correctRT_using_pg=True, verbose=True)
        alignment.alignBestCluster_legacy(self.multipeptides, tree_mapped, self.tr_data)

        # We should have 4 peakgroups
        prec1 = self.mpep
        self.assertEqual(len(prec1.get_selected_peakgroups()), 4)

        # Check that we have all the correct ones (1,2,4,8)
        self.assertEqual(set(['peakgroup8', 'peakgroup2', 'peakgroup4', 'peakgroup1']), 
                         set([p.get_feature_id() for p in prec1.get_selected_peakgroups()]))

    def test_alignBestCluster_1(self):
        """Test the best cluster align

        This is now using no correction of the alignment thread by using the
        found peakgroup. In this case it means that after finding the second
        peakgroup at 112 s, the search RT for run 2 is still at 106 seconds
        which gets mapped to 112 seconds in run 3 (but the next pg is at 120s,
        too far for 7 seconds tolerance).
        """

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(algo.getDistanceMatrix(self.exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.exp.runs[a].get_id(), self.exp.runs[b].get_id()) for a,b in tree]

        alignment = algo.TreeConsensusAlignment(max_rt_diff = 6, fdr_cutoff = 0.1, aligned_fdr_cutoff = 0.25, correctRT_using_pg=False)
        alignment.alignBestCluster_legacy(self.multipeptides, tree_mapped, self.tr_data)

        # Now only 2 peakgroups should be selected
        prec1 = self.mpep
        self.assertEqual(len(prec1.get_selected_peakgroups()), 2)

        # Check that we have all the correct ones (only 1,2)
        self.assertEqual(set(['peakgroup2', 'peakgroup1']), 
                         set([p.get_feature_id() for p in prec1.get_selected_peakgroups()]))

    def test_alignBestCluster_2(self):
        """Test the best cluster align

        This is now using no correction of the alignment thread by using the
        found peakgroup (e.g. no correction of the threading). 

          - Run1 : 100s     [threadRT = 100s] 
          - Run2 : 112s     [threadRT = 106s]
          - Run3 : 120s     [threadRT = 112s]
          - Run4 : xxx      [threadRT = 118s]
          - Run5 : 139s     [threadRT = 124s]

        By using a larger tolerance of 15s, we can still manage to find all the correct peakgroups
        """

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(algo.getDistanceMatrix(self.exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.exp.runs[a].get_id(), self.exp.runs[b].get_id()) for a,b in tree]

        alignment = algo.TreeConsensusAlignment(max_rt_diff = 15, fdr_cutoff = 0.1, aligned_fdr_cutoff = 0.25, correctRT_using_pg=False)
        alignment.alignBestCluster_legacy(self.multipeptides, tree_mapped, self.tr_data)

        # Now only 2 peakgroups should be selected
        prec1 = self.mpep
        self.assertEqual(len(prec1.get_selected_peakgroups()), 4)

        # Check that we have all the correct ones (1,2,4,8)
        self.assertEqual(set(['peakgroup8', 'peakgroup2', 'peakgroup4', 'peakgroup1']), 
                         set([p.get_feature_id() for p in prec1.get_selected_peakgroups()]))

    def test_alignAllCluster_1(self):
        """Test the best cluster align
        
        This is using the best possible conditions with only 7 seconds retention time cutoff

          - Run1 : 100s     [threadRT = 100s] 
          - Run2 : 112s     [threadRT = 106s]
          - Run3 : 120s     [threadRT = 118s]
          - Run4 : xxx      [threadRT = 126s]  (should be around 130s)
          - Run5 : 139s     [threadRT = 133s]
        """

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(algo.getDistanceMatrix(self.exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.exp.runs[a].get_id(), self.exp.runs[b].get_id()) for a,b in tree]

        alignment = algo.TreeConsensusAlignment(max_rt_diff = 6, fdr_cutoff = 0.1, aligned_fdr_cutoff = 0.25, correctRT_using_pg=True, verbose=True)
        alignment.alignAllCluster(self.multipeptides, tree_mapped, self.tr_data)

        # We should have 4 peakgroups
        prec1 = self.mpep
        self.assertEqual(len(prec1.get_selected_peakgroups()), 4)

        # Check that we have all the correct ones (1,2,4,8)
        self.assertEqual(set(['peakgroup8', 'peakgroup2', 'peakgroup4', 'peakgroup1']), 
                         set([p.get_feature_id() for p in prec1.get_selected_peakgroups()]))

    def test_alignAllCluster_2(self):
        """Test the best cluster align
        
        This is using the best possible conditions with only 7 seconds retention time cutoff

        Cluster 1:
          - Run1 : 100s     [threadRT = 100s] 
          - Run2 : 112s     [threadRT = 106s]
          - Run3 : 120s     [threadRT = 118s]
          - Run4 : xxx      [threadRT = 126s]  (should be around 130s)
          - Run5 : 139s     [threadRT = 133s]


        Cluster 2:
          - Run3 : 130s     [threadRT = 130s]
          - Run4 : 145s     [threadRT = 139s]
          - Run2 : 130s     [threadRT = 122s]

        """

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(algo.getDistanceMatrix(self.exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.exp.runs[a].get_id(), self.exp.runs[b].get_id()) for a,b in tree]

        alignment = algo.TreeConsensusAlignment(max_rt_diff = 9, fdr_cutoff = 0.2, aligned_fdr_cutoff = 0.25, correctRT_using_pg=True, verbose=True)
        alignment.alignAllCluster(self.multipeptides, tree_mapped, self.tr_data)

        # We should have 4 peakgroups selected and 7 peakgroups in clusters
        prec1 = self.mpep
        self.assertEqual(len(prec1.get_selected_peakgroups()), 4)
        self.assertEqual(len([pg for pep in prec1.getAllPeptides() for pg in pep.get_all_peakgroups()]), 7)

        # Check that we have all the correct ones (1,2,4,8)
        self.assertEqual(set(['peakgroup8', 'peakgroup2', 'peakgroup4', 'peakgroup1']), 
                         set([p.get_feature_id() for p in prec1.get_selected_peakgroups()]))

        pg_cluster1 = [pg for pep in prec1.getAllPeptides() for pg in pep.get_all_peakgroups() if pg.get_cluster_id() == 1]
        pg_cluster2 = [pg for pep in prec1.getAllPeptides() for pg in pep.get_all_peakgroups() if pg.get_cluster_id() == 2]

        # Check the two individual clusters
        self.assertEqual(len(pg_cluster1), 4)
        self.assertEqual(len(pg_cluster2), 3)
        self.assertEqual(set(['peakgroup8', 'peakgroup2', 'peakgroup4', 'peakgroup1']), 
                         set([p.get_feature_id() for p in pg_cluster1]))
        self.assertEqual(set(['peakgroup7', 'peakgroup3', 'peakgroup5']), 
                         set([p.get_feature_id() for p in pg_cluster2]))

if __name__ == '__main__':
    unittest.main()
