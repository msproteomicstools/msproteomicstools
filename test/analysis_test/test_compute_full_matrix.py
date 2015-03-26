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

from __future__ import print_function
import unittest
import subprocess as sub
import os
from nose.plugins.attrib import attr

class TestComputeFullMatrix(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "analysis")

    def exact_diff(self, name1, name2, sep="\t", header_exclude = []):
        """ Check whether two (csv/tsv) files are almost equal, allowing for numerical inaccuracies only."""

        f1 = open(name1, "r")
        f2 = open(name2, "r")
        header_nr_exclude = []
        for i,(l1,l2) in enumerate(zip(f1,f2)):

            # Iterate through all fields (split by seperator)
            for j,(field1,field2) in enumerate(zip(l1.split(sep),l2.split(sep))):

                # Exclude certain requested columns
                if i == 0 and field1 in header_exclude:
                    header_nr_exclude.append(j)
                elif j in header_nr_exclude:
                    continue

                # Check all other fields for equality
                try:
                    self.assertAlmostEqual(float(field1),float(field2) )
                except ValueError:
                    self.assertEqual(field1,field2)

    def test_1_full_matrix(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "compute_full_matrix.py")
        filename = os.path.join(self.datadir, "computeFullMatrix_1_input.csv")
        expected_matrix_outcome = os.path.join(self.datadir, "computeFullMatrix_1_output.tsv")
        tmpfilename = "fullMatrix_1.out.tmp.tsv"

        args = "--in %s --out_matrix %s --aligner_mscore_threshold 1.0 --output_method full" % (filename, tmpfilename)
        cmd = "python %s %s" % (script, args)
        print(cmd)
        sub.check_output(cmd,shell=True)
        
        self.exact_diff(tmpfilename, expected_matrix_outcome)

        os.remove(tmpfilename)

if __name__ == '__main__':
    unittest.main()
