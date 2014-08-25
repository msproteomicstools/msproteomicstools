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
import subprocess as sub
import os
from nose.plugins.attrib import attr

from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader, inferMapping
from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment
from msproteomicstoolslib.algorithms.alignment.AlignmentMST import getDistanceMatrix
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree
import msproteomicstoolslib.math.Smoothing as smoothing
import msproteomicstoolslib.format.TransformationCollection as transformations
from msproteomicstoolslib.algorithms.alignment.BorderIntegration import \
        integrationBorderShortestPath, integrationBorderShortestDistance, integrationBorderReference
import msproteomicstoolslib.algorithms.alignment.AlignmentHelper as helper
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner

class TestAlignment(unittest.TestCase):

    def setUp(self):
        # Set up dirs
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "analysis")

        # Set up files
        peakgroups_file = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        mzml_file = os.path.join(self.datadir, "imputeValues/r004_small/split_olgas_otherfile.chrom.mzML")

        # Parameters
        self.initial_alignment_cutoff = 0.0001
        fdr_cutoff_all_pg = 1.0
        max_rt_diff = 30

        # Read input
        reader = SWATHScoringReader.newReader([peakgroups_file], "openswath", readmethod="complete")
        self.new_exp = MRExperiment()
        self.new_exp.runs = reader.parse_files()
        self.multipeptides = self.new_exp.get_all_multipeptides(fdr_cutoff_all_pg, verbose=False)

        # Align all against all
        self.tr_data = transformations.LightTransformationData()
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        for run_0 in self.new_exp.runs:
            for run_1 in self.new_exp.runs:
                helper.addDataToTrafo(self.tr_data, run_0, run_1, spl_aligner, self.multipeptides, "linear", 30)

        # Select two interesting peptides
        pepname = "21517_C[160]NVVISGGTGSGK/2_run0 0 0"
        self.current_mpep1 = [m for m in self.multipeptides if m.get_peptides()[0].get_id() == pepname][0]

        pepname = "26471_GYEDPPAALFR/2_run0 0 0"
        self.current_mpep2 = [m for m in self.multipeptides if m.get_peptides()[0].get_id() == pepname][0]

    def test_shortestDistance_1(self):

        rid = "0_0"
        dist_matrix = getDistanceMatrix(self.new_exp, self.multipeptides, self.initial_alignment_cutoff)

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.get_peptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        rmap = dict([(r.get_id(),i) for i,r in enumerate(self.new_exp.runs) ])
        border_l, border_r = integrationBorderShortestDistance(selected_pg, 
            rid, self.tr_data, dist_matrix, rmap)

        # Direct transformation from 0_2 to 0_0
        self.assertAlmostEqual(border_l, self.tr_data.getTrafo("0_2", "0_0").predict([ 240.0 ])[0])
        self.assertAlmostEqual(border_r, self.tr_data.getTrafo("0_2", "0_0").predict([ 260.0 ])[0])

        self.assertAlmostEqual(border_l, 77.992277992277934)
        self.assertAlmostEqual(border_r, 84.1698841699)

    def test_shortestPath_1(self):

        rid = "0_0"
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, self.initial_alignment_cutoff))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.get_peptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderShortestPath(selected_pg, 
            rid, self.tr_data, tree_mapped)

        # Direct transformation from 0_2 to 0_0
        self.assertAlmostEqual(border_l, self.tr_data.getTrafo("0_2", "0_0").predict([ 240.0 ])[0])
        self.assertAlmostEqual(border_r, self.tr_data.getTrafo("0_2", "0_0").predict([ 260.0 ])[0])

        self.assertAlmostEqual(border_l, 77.992277992277934)
        self.assertAlmostEqual(border_r, 84.1698841699)

    def test_reference_1(self):

        rid = "0_0"
        self.tr_data.reference = "0_2" # set reference run to 0_2

        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, self.initial_alignment_cutoff))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.get_peptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderReference(self.new_exp, selected_pg, 
            rid, self.tr_data, "median")

        # Direct transformation from 0_2 to 0_0
        self.assertAlmostEqual(border_l, self.tr_data.getTrafo("0_2", "0_0").predict([ 240.0 ])[0])
        self.assertAlmostEqual(border_r, self.tr_data.getTrafo("0_2", "0_0").predict([ 260.0 ])[0])

        self.assertAlmostEqual(border_l, 77.992277992277934)
        self.assertAlmostEqual(border_r, 84.1698841699)

        border_l, border_r = integrationBorderReference(self.new_exp, selected_pg, 
            rid, self.tr_data, "mean")
        self.assertAlmostEqual(border_l, 77.992277992277934)
        self.assertAlmostEqual(border_r, 84.1698841699)
        border_l, border_r = integrationBorderReference(self.new_exp, selected_pg, 
            rid, self.tr_data, "max_width")
        self.assertAlmostEqual(border_l, 77.992277992277934)
        self.assertAlmostEqual(border_r, 84.1698841699)

        self.assertRaises(Exception, integrationBorderReference, self.new_exp, selected_pg, 
            rid, self.tr_data, "dummy")

    def test_shortestDistance_2(self):

        rid = "0_1"
        dist_matrix = getDistanceMatrix(self.new_exp, self.multipeptides, self.initial_alignment_cutoff)

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.get_peptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        rmap = dict([(r.get_id(),i) for i,r in enumerate(self.new_exp.runs) ])
        border_l, border_r = integrationBorderShortestDistance(selected_pg, 
            rid, self.tr_data, dist_matrix, rmap)

        # Shortest distance means that we transformed directly from 0_2 to 0_1
        self.assertAlmostEqual(border_l, self.tr_data.getTrafo("0_2", "0_1").predict([ 240.0 ])[0])
        self.assertAlmostEqual(border_r, self.tr_data.getTrafo("0_2", "0_1").predict([ 260.0 ])[0])

        self.assertAlmostEqual(border_l, 168.03088803088787)
        self.assertAlmostEqual(border_r, 183.32046332)

    def test_shortestPath_2(self):

        rid = "0_1"
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, self.initial_alignment_cutoff))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.get_peptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderShortestPath(selected_pg, 
            rid, self.tr_data, tree_mapped)

        # Shortest path means that we transformed from 0_2 to 0_0 and then to 0_1
        self.assertAlmostEqual(border_l, 
                               self.tr_data.getTrafo("0_0", "0_1").predict(
                                 self.tr_data.getTrafo("0_2", "0_0").predict([ 240.0 ]) 
                               ))
        self.assertAlmostEqual(border_r, 
                               self.tr_data.getTrafo("0_0", "0_1").predict(
                                 self.tr_data.getTrafo("0_2", "0_0").predict([ 260.0 ]) 
                               ))

        self.assertAlmostEqual(border_l, 187.18146718146681)
        self.assertAlmostEqual(border_r, 202.00772200772167)

    def test_reference_2(self):

        rid = "0_1"
        self.tr_data.reference = "0_0" # set reference run to 0_0

        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, self.initial_alignment_cutoff))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.get_peptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderReference(self.new_exp, selected_pg, 
            rid, self.tr_data, "median")

        # Shortest path means that we transformed from 0_2 to 0_0 and then to 0_1
        self.assertAlmostEqual(border_l, 
                               self.tr_data.getTrafo("0_0", "0_1").predict(
                                 self.tr_data.getTrafo("0_2", "0_0").predict([ 240.0 ]) 
                               ))
        self.assertAlmostEqual(border_r, 
                               self.tr_data.getTrafo("0_0", "0_1").predict(
                                 self.tr_data.getTrafo("0_2", "0_0").predict([ 260.0 ]) 
                               ))

        self.assertAlmostEqual(border_l, 187.18146718146681)
        self.assertAlmostEqual(border_r, 202.00772200772167)

if __name__ == '__main__':
    unittest.main()
