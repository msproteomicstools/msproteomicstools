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

from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader
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
        self.current_mpep1 = [m for m in self.multipeptides if m.getAllPeptides()[0].get_id() == pepname][0]

        pepname = "26471_GYEDPPAALFR/2_run0 0 0"
        self.current_mpep2 = [m for m in self.multipeptides if m.getAllPeptides()[0].get_id() == pepname][0]

    def test_shortestDistance_1(self):

        rid = "0_0"
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        dist_matrix = getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner)

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
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
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
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

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
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
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        dist_matrix = getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner)

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
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
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderShortestPath(selected_pg, 
            rid, self.tr_data, tree_mapped)

        # Shortest path means that we transformed from 0_2 to 0_1
        self.assertAlmostEqual(border_l, self.tr_data.getTrafo("0_2", "0_1").predict( [ 240.0 ] ))
        self.assertAlmostEqual(border_r, self.tr_data.getTrafo("0_2", "0_1").predict( [ 260.0 ] ))

        self.assertAlmostEqual(border_l, 168.03088803088787)
        self.assertAlmostEqual(border_r, 183.32046332046318)

    def test_shortestPath_3(self):

        rid = "0_1"
        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep2.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderShortestPath(selected_pg, 
            rid, self.tr_data, tree_mapped)

        # Shortest path means that we transformed from 0_0 to 0_2 and then to 0_1
        self.assertAlmostEqual(border_l, 
                               self.tr_data.getTrafo("0_2", "0_1").predict(
                                 self.tr_data.getTrafo("0_0", "0_2").predict([ 600 ]) 
                               ))
        self.assertAlmostEqual(border_r, 
                               self.tr_data.getTrafo("0_2", "0_1").predict(
                                 self.tr_data.getTrafo("0_0", "0_2").predict([ 700 ]) 
                               ))

        self.assertAlmostEqual(border_l, 1452.355212355212)
        self.assertAlmostEqual(border_r, 1696.9884169884167)

    def test_reference_2(self):

        rid = "0_1"
        self.tr_data.reference = "0_0" # set reference run to 0_0

        spl_aligner = SplineAligner(self.initial_alignment_cutoff)
        tree = MinimumSpanningTree(getDistanceMatrix(self.new_exp, self.multipeptides, spl_aligner))
        tree_mapped = [(self.new_exp.runs[a].get_id(), self.new_exp.runs[b].get_id()) for a,b in tree]

        # Select peakgroups, compute left/right border
        selected_pg = [pg for p in self.current_mpep1.getAllPeptides() for pg in p.get_all_peakgroups() if pg.get_cluster_id() == 1]
        border_l, border_r = integrationBorderReference(self.new_exp, selected_pg, 
            rid, self.tr_data, "median")

        # Reference 0_0 means that we transformed from 0_2 to 0_0 and then to 0_1
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

class MatrixOutputWriters(unittest.TestCase):

    def setUp(self):

        # Set up dirs
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "analysis")

        # Set up files
        peakgroups_file = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        fdr_cutoff_all_pg = 1.0

        # Read input
        reader = SWATHScoringReader.newReader([peakgroups_file], "openswath", readmethod="complete")
        self.exp = MRExperiment()
        self.exp.runs = reader.parse_files()
        self.multipeptides = self.exp.get_all_multipeptides(fdr_cutoff_all_pg, verbose=False)

        # Set up files nr2 
        peakgroups_file = os.path.join(self.datadir, "feature_alignment_7_openswath_input.csv")
        reader = SWATHScoringReader.newReader([peakgroups_file], "openswath", readmethod="complete")
        self.exp2 = MRExperiment()
        self.exp2.runs = reader.parse_files()
        self.multipeptides2 = self.exp2.get_all_multipeptides(fdr_cutoff_all_pg, verbose=False)

        # Select the best peakgroup per peptide and select it for writing out
        fdr_cutoff = 0.01
        for mpep in self.multipeptides2:
            for prgr in mpep.getAllPeptides():
                minpg = min( [(pg.get_fdr_score(), pg) for pg in prgr.peakgroups] )
                if minpg[0] < fdr_cutoff:
                    minpg[1].select_this_peakgroup()

    def test_matrix_out_1(self):
        """Test the output matrix writers"""

        import msproteomicstoolslib.algorithms.alignment.AlignmentHelper as helper

        tmpfile = "tmp.output.csv"
        helper.write_out_matrix_file(tmpfile, self.exp.runs, self.multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

        tmpfile = "tmp.output.tsv"
        helper.write_out_matrix_file(tmpfile, self.exp.runs, self.multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

        tmpfile = "tmp.output.xls"
        helper.write_out_matrix_file(tmpfile, self.exp.runs, self.multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)
        tmpfile = "tmp.output.xlsx"
        helper.write_out_matrix_file(tmpfile, self.exp.runs, self.multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

    def test_matrix_out_2(self):
        """Test the output matrix writers"""

        import msproteomicstoolslib.algorithms.alignment.AlignmentHelper as helper

        runs = self.exp2.runs
        multipeptides = self.multipeptides2

        tmpfile = "tmp.output.csv"
        helper.write_out_matrix_file(tmpfile, runs, multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

        tmpfile = "tmp.output.tsv"
        helper.write_out_matrix_file(tmpfile, runs, multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

        tmpfile = "tmp.output.xls"
        helper.write_out_matrix_file(tmpfile, runs, multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)
        tmpfile = "tmp.output.xlsx"
        helper.write_out_matrix_file(tmpfile, runs, multipeptides, 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

    def test_matrix_out_3(self):
        """Test the output matrix writers"""

        import msproteomicstoolslib.algorithms.alignment.AlignmentHelper as helper

        runs = self.exp2.runs
        multipeptides = self.multipeptides2[0]  

        # Try to mess up the assumption of one selected peakgroup per run
        pep1 = multipeptides.getAllPeptides()[0]
        pep2 = multipeptides.getAllPeptides()[0]
        for pg in pep1.peakgroups:
            pg.select_this_peakgroup()
        for pg in pep2.peakgroups:
            pg.peptide = pep1
            pg.select_this_peakgroup()

        tmpfile = "tmp.output.csv"
        self.assertRaises(Exception, helper.write_out_matrix_file, tmpfile, runs, [ multipeptides ] , 0.0, 
                                     style="full", write_requant=False)
        os.remove(tmpfile)

if __name__ == '__main__':
    unittest.main()

