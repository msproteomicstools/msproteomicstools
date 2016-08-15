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

class TestNoiseIntegration(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "analysis")

    def exact_diff(self, name1, name2):
        """ Check whether two (csv/tsv) files are almost equal, allowing for numerical inaccuracies only."""
        f1 = open(name1, "r")
        f2 = open(name2, "r")
        try:
            for l1,l2 in zip(f1,f2):
                for field1,field2 in zip(l1.split(),l2.split()):
                    try:
                        self.assertAlmostEqual(float(field1),float(field2) )
                    except ValueError:
                        self.assertEqual(field1,field2)
        except AssertionError as e:
            print("Found non-identical files, for investigation please use:")
            print("  diff %s %s" % (name1, name2))
            raise e

    @attr('slow')
    def test_1_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_1_input.csv")
        tr_f1 = os.path.join(self.datadir, "imputeValues/r003_small/transformation-0_0-0_0.tr")
        tr_f2 = os.path.join(self.datadir, "imputeValues/r004_small/transformation-0_1-0_0.tr")
        expected_outcome = os.path.join(self.datadir, "imputeValues_1_.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_1_output_matrix.csv")
        tmpfilename = "imputeValues_1.out.tmp"
        tmpfilename_matrix = "imputeValues_1.out.tmp_matrix.tsv"

        args = "--in %s %s --peakgroups_infile %s --out %s --out_matrix %s --border_option median --matrix_output_method RT --test" % (
            tr_f1, tr_f2, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        # self.exact_diff(tmpfilename_ids, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    @attr('slow')
    def test_2_cache_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_1_input.csv")
        tr_f1 = os.path.join(self.datadir, "imputeValues/r003_small/transformation-0_0-0_0.tr")
        tr_f2 = os.path.join(self.datadir, "imputeValues/r004_small/transformation-0_1-0_0.tr")
        expected_outcome = os.path.join(self.datadir, "imputeValues_1_.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_1_output_matrix.csv")
        tmpfilename = "imputeValues_2.out.tmp"
        tmpfilename_matrix = "imputeValues_2.out.tmp_matrix.tsv"

        # We should get the same results if we cache the chromatograms in memory
        args = "--in %s %s --peakgroups_infile %s --out %s --out_matrix %s --cache_in_memory --border_option median --matrix_output_method RT --test" % (
            tr_f1, tr_f2, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        # self.exact_diff(tmpfilename_ids, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    def test_3_singleShortestPath_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_1_input.csv")
        mzml_file = os.path.join(self.datadir, "imputeValues/r004_small/split_olgas_otherfile.chrom.mzML")
        expected_outcome = os.path.join(self.datadir, "imputeValues_3_output.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_3_output_matrix.csv")
        tmpfilename = "imputeValues_3.out.tmp"
        tmpfilename_matrix = "imputeValues_3.out.tmp_matrix.tsv"

        args = "--do_single_run %s --peakgroups_infile %s  --out %s --out_matrix %s\
                --realign_runs None --method singleShortestPath --matrix_output_method RT --test" % (
            mzml_file, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    def test_4_singleClosestRun_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_1_input.csv")
        mzml_file = os.path.join(self.datadir, "imputeValues/r004_small/split_olgas_otherfile.chrom.mzML")
        expected_outcome = os.path.join(self.datadir, "imputeValues_3_output.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_3_output_matrix.csv")
        tmpfilename = "imputeValues_4.out.tmp"
        tmpfilename_matrix = "imputeValues_4.out.tmp_matrix.tsv"

        args = "--do_single_run %s --peakgroups_infile %s  --out %s --out_matrix %s\
                --realign_runs None --method singleClosestRun --matrix_output_method RT --test" % (
            mzml_file, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    def test_5_singleClosestRun_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        mzml_file = os.path.join(self.datadir, "imputeValues/r004_small/split_olgas_otherfile.chrom.mzML")
        expected_outcome = os.path.join(self.datadir, "imputeValues_5_output.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_5_output_matrix.csv")
        tmpfilename = "imputeValues_5.out.tmp"
        tmpfilename_matrix = "imputeValues_5.out.tmp_matrix.tsv"

        args = "--do_single_run %s --peakgroups_infile %s  --out %s --out_matrix %s\
                --realign_runs linear --method singleClosestRun --matrix_output_method RT --test" % (
            mzml_file, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    def test_6_singleShortestPath_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_5_input.csv")
        mzml_file = os.path.join(self.datadir, "imputeValues/r004_small/split_olgas_otherfile.chrom.mzML")
        expected_outcome = os.path.join(self.datadir, "imputeValues_6_output.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_6_output_matrix.csv")
        tmpfilename = "imputeValues_6.out.tmp"
        tmpfilename_matrix = "imputeValues_6.out.tmp_matrix.tsv"

        args = "--do_single_run %s --peakgroups_infile %s  --out %s --out_matrix %s\
                --realign_runs linear --method singleShortestPath --matrix_output_method RT --test" % (
            mzml_file, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    @attr('slow')
    def test_7_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_7_input.csv")
        tr_f1 = os.path.join(self.datadir, "imputeValues/r003_small/transformation-0_0-0_0.tr")
        tr_f2 = os.path.join(self.datadir, "imputeValues/r004_small/transformation-0_1-0_0.tr")
        expected_outcome = os.path.join(self.datadir, "imputeValues_7_output.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_7_output_matrix.csv")
        tmpfilename = "imputeValues_7.out.tmp"
        tmpfilename_matrix = "imputeValues_7.out.tmp_matrix.tsv"

        args = "--in %s %s --peakgroups_infile %s --out %s --out_matrix %s --border_option median --matrix_output_method RT --disable_isotopic_transfer --test" % (
            tr_f1, tr_f2, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

    @attr('slow')
    def test_8_requantAlignedValues(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "requantAlignedValues.py")
        filename = os.path.join(self.datadir, "imputeValues/imputeValues_8_input.csv")
        tr_f1 = os.path.join(self.datadir, "imputeValues/r003_small/transformation-0_0-0_0.tr")
        tr_f2 = os.path.join(self.datadir, "imputeValues/r004_small/transformation-0_1-0_0.tr")
        expected_outcome = os.path.join(self.datadir, "imputeValues_8_output.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "imputeValues_8_output_matrix.csv")
        tmpfilename = "imputeValues_8.out.tmp"
        tmpfilename_matrix = "imputeValues_8.out.tmp_matrix.tsv"

        args = "--in %s %s --peakgroups_infile %s --out %s --out_matrix %s --border_option median --matrix_output_method RT --test" % (
            tr_f1, tr_f2, filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename_matrix)

if __name__ == '__main__':
    unittest.main()
