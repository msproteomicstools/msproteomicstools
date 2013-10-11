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

class TestFeatureAlignment(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "analysis")

    def exact_diff(self, name1, name2):
        """ Check whether two (csv/tsv) files are almost equal, allowing for numerical inaccuracies only."""
        f1 = open(name1, "r")
        f2 = open(name2, "r")
        for l1,l2 in zip(f1,f2):
            for field1,field2 in zip(l1.split(),l2.split()):
                try:
                    self.assertAlmostEqual(float(field1),float(field2) )
                except ValueError:
                    self.assertEqual(field1,field2)


    def test_1_featureAlignment_openswath(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "feature_alignment.py")
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        expected_outcome = os.path.join(self.datadir, "feature_alignment_1_openswath_output_cluster_ids.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "feature_alignment_1_openswath_output_matrix.csv")
        tmpfilename = "featureAlignment_1.out.tmp"
        tmpfilename_ids = "featureAlignment_1.out.tmp_idsonly.csv"
        tmpfilename_matrix = "featureAlignment_1.out.tmp_matrix.csv"

        args = "--in %s --out %s --out_ids %s --out_matrix %s --method best_cluster_score --max_fdr_quality 0.4" % (filename, tmpfilename, tmpfilename_ids, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename_ids, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename)
        os.remove(tmpfilename_ids)
        os.remove(tmpfilename_matrix)

    def test_2_featureAlignment_openswath_best_overall(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "feature_alignment.py")
        filename = os.path.join(self.datadir, "feature_alignment_openswath_input_1.csv")
        expected_outcome = os.path.join(self.datadir, "feature_alignment_2_output_1_ids.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "feature_alignment_2_output_2_matrix.csv")
        tmpfilename = "featureAlignment_2.out.tmp"
        tmpfilename_ids = "featureAlignment_2.out.tmp_idsonly.csv"
        tmpfilename_matrix = "featureAlignment_2.out.tmp_matrix.csv"

        args = "--in %s --out %s --out_ids %s --out_matrix %s --method best_overall --max_fdr_quality 0.4" % (filename, tmpfilename, tmpfilename_ids, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename_ids, expected_outcome)
        self.exact_diff(tmpfilename_matrix, expected_matrix_outcome)

        os.remove(tmpfilename)
        os.remove(tmpfilename_ids)
        os.remove(tmpfilename_matrix)

    def test_3_featureAlignment_openswath_alignment(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "feature_alignment.py")
        filename = os.path.join(self.datadir, "feature_alignment_3_openswath_input.csv")
        expected_outcome = os.path.join(self.datadir, "feature_alignment_3_openswath_output_cluster_ids.csv")
        tmpfilename = "featureAlignment_3.out.tmp"
        tmpfilename_ids = "featureAlignment_3.out.tmp_idsonly.csv"

        args = "--in %s --out %s --out_ids %s --realign_runs --method best_cluster_score --max_fdr_quality 0.4" % (filename, tmpfilename, tmpfilename_ids)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename_ids, expected_outcome)

        os.remove(tmpfilename)
        os.remove(tmpfilename_ids)

    def test_4_featureAlignment_openswath_alignment_scikit(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "feature_alignment.py")
        filename = os.path.join(self.datadir, "feature_alignment_3_openswath_input.csv")
        expected_outcome = os.path.join(self.datadir, "feature_alignment_3_openswath_output_cluster_ids.csv")
        tmpfilename = "featureAlignment_4.out.tmp"
        tmpfilename_ids = "featureAlignment_4.out.tmp_idsonly.csv"

        args = "--in %s --out %s --out_ids %s --realign_runs --use_scikit --method best_cluster_score --max_fdr_quality 0.4" % (filename, tmpfilename, tmpfilename_ids)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename_ids, expected_outcome)

        os.remove(tmpfilename)
        os.remove(tmpfilename_ids)

    def test_5_featureAlignment_peakview(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "feature_alignment.py")
        filename = os.path.join(self.datadir, "feature_alignment_peakview_input_2.csv")
        expected_outcome = os.path.join(self.datadir, "feature_alignment_5_peakview_output_matrix.csv")
        tmpfilename = "featureAlignment_5.out.tmp"
        tmpfilename_matrix = "featureAlignment_5.out.tmp_matrix.csv"

        args = "--in %s --out %s --out_matrix %s --file_format peakview  --outlier_thresh 5 --max_fdr_quality 0.0001 --fdr_cutoff 0.000000001 --method best_cluster_score" % (filename, tmpfilename, tmpfilename_matrix)
        cmd = "python %s %s" % (script, args)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename_matrix, expected_outcome)

        os.remove(tmpfilename_matrix)

if __name__ == '__main__':
    unittest.main()
