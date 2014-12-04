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
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""

import unittest
import os

import msproteomicstoolslib.format.SWATHScoringReader as reader

def doTest(self, runs):
    self.assertEqual(len(runs), 4)
    self.assertEqual(runs[0].get_id(), "0_0")
    self.assertEqual(runs[1].get_id(), "1_0")
    self.assertEqual(runs[2].get_id(), "2_0")
    self.assertEqual(runs[3].get_id(), "3_0")

    self.assertEqual(runs[0].get_openswath_filename(), "split_hroest_K120808_combined.featureXML")
    self.assertIsNone(runs[0].get_aligned_filename())

    self.assertEqual(len(list(runs[0])), 1)
    self.assertEqual(len(list(runs[1])), 1)
    self.assertEqual(len(list(runs[2])), 1)
    self.assertEqual(len(list(runs[3])), 1)

    self.assertIsNotNone(runs[0].getPrecursorGroup("17365_TLDTAAEKIVETATR/3_run0"))
    self.assertIsNone(runs[0].getPrecursorGroup("19365_TLDTAAEKIVETATR/3_run0"))

class TestUnitScoringReaderOpenSWATH(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")

    def test_newReader(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "openswath", "complete")
        self.assertTrue(True)

    def test_parse_files_osw_complete(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "openswath", "complete")
        runs = r.parse_files(False)

        doTest(self, runs)

    def test_parse_files_osw_min(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "openswath", "minimal")
        runs = r.parse_files(False)

        doTest(self, runs)

    def test_parse_files_osw_min_exp(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "openswath", "minimal")
        runs = r.parse_files(True)

        doTest(self, runs)

    def test_parse_files_osw_gui(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "openswath", "gui")
        runs = r.parse_files(False)

        doTest(self, runs)

    def test_parse_files_clusterId(self):
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        r = reader.SWATHScoringReader.newReader([filename], "openswath", "minimal")
        runs = r.parse_files(False)

        prgroup = runs[2].getPrecursorGroup("21517_C[160]NVVISGGTGSGK/2_run0 0 0")
        # The precursor is now the one and only precursor present
        self.assertEqual(len(prgroup.getAllPrecursors()), 1)
        p = prgroup.getAllPrecursors()[0]

        all_pg = list(p.getAllPeakgroups())
        cl_pg = list(p.getClusteredPeakgroups())

        self.assertEqual(len(cl_pg), 2)
        self.assertEqual(len(all_pg), 2)
        self.assertEqual(all_pg[0].get_cluster_id(), 1)
        self.assertEqual(all_pg[1].get_cluster_id(), 2)

class TestUnitScoringReaderPeakView(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")

    def test_newReader(self):
        filename = os.path.join(self.datadir, "feature_alignment_peakview_input_2.csv")
        r = reader.SWATHScoringReader.newReader([filename], "peakview", "minimal")
        self.assertTrue(True)

    def test_parse_files_peakview_min(self):
        filename = os.path.join(self.datadir, "feature_alignment_peakview_input_2.csv")
        r = reader.SWATHScoringReader.newReader([filename], "peakview", "minimal")
        runs = r.parse_files(False)

        self.assertEqual(len(runs), 5)
        self.assertEqual(runs[0].get_id(), 'ywu_L121218_002_SW.wiff (sample 1)_0')

        self.assertIsNone(runs[0].get_openswath_filename())
        self.assertIsNone(runs[0].get_aligned_filename())

        self.assertEqual(len(list(runs[0])), 2)
        self.assertEqual(len(list(runs[1])), 2)
        self.assertEqual(len(list(runs[2])), 2)
        self.assertEqual(len(list(runs[3])), 2)
        self.assertEqual(len(list(runs[4])), 2)

        self.assertIsNotNone(runs[0].getPrecursorGroup('LIGNMALLPLR+2'))
        self.assertIsNotNone(runs[0].getPrecursorGroup('LIGNMALLPLR+3'))
        self.assertIsNone(runs[0].getPrecursorGroup('LIGNMALLPLR+4'))

class TestUnitScoringReadermProphet(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")

    def test_newReader(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "mprophet", "complete")
        self.assertTrue(True)

    def test_parse_files_mprophet_complete(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "mprophet", "complete")
        runs = r.parse_files(False)

        doTest(self, runs)

    def test_parse_files_mprophet_min(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        r = reader.SWATHScoringReader.newReader([filename], "mprophet", "minimal")
        runs = r.parse_files(False)

        doTest(self, runs)

class TestUnitScoring(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")

    def test_inferMapping(self):
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        rawfiles = ["split_olgas_K121102_003_SW_PA1_earlyexp", "somethign_else", "split_olgas_otherfile"]
        mapping = {}
        precursors_mapping = {}
        sequences_mapping = {}
        protein_mapping = {}
        reader.inferMapping(rawfiles, [filename], mapping, verbose=True, throwOnMismatch=False, 
            precursors_mapping = precursors_mapping, sequences_mapping = sequences_mapping, protein_mapping = protein_mapping)
        expected = {'0_2': ['somethign_else'], '0_0': ['split_olgas_K121102_003_SW_PA1_earlyexp'], '0_1': ['split_olgas_otherfile']}
        self.assertEqual(mapping, expected)

    def test_inferThrow(self):
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        rawfiles = ["split_olgas_K121102_003_SW_PA1_earlyexp", "split_olgas_otherfile"]
        mapping = {}
        precursors_mapping = {}
        sequences_mapping = {}
        protein_mapping = {}
        self.assertRaises(Exception, reader.inferMapping, rawfiles, [filename], mapping, precursors_mapping, sequences_mapping, protein_mapping, False, True)
        expected = {'0_0': ['split_olgas_K121102_003_SW_PA1_earlyexp'], '0_1': ['split_olgas_otherfile']}
        self.assertEqual(mapping, expected)

    def test_inferWrongFile(self):
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        self.assertRaises(Exception, reader.inferMapping, [], [filename], {}, False, True)

class TestUnitSWATHScoringReader(unittest.TestCase):

    def setUp(self):
        pass

    def test(self):
        self.assertRaises(Exception, reader.SWATHScoringReader)
        self.assertRaises(Exception, reader.SWATHScoringReader.newReader, ["test"], "DoesntExist", "complete")

if __name__ == '__main__':
    unittest.main()

