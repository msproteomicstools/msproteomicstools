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

from msproteomicstoolslib.data_structures.Run import Run

class MockPrecursorGroup():

    def __init__(self, id_):
        self.id_ = id_

    def get_best_peakgroup(self):
        return "42"

    # Dummy function - instead of yielding a second level Peptide, it simply yields itself
    def __iter__(self):
        yield self

class TestUnitRun(unittest.TestCase):

    def setUp(self):
        pass

    def test_createRun(self):
        r = Run([], {}, "run1", "file1.txt", filename="file1.csv", aligned_filename="file1.tsv")
        self.assertTrue(True)
        self.assertEqual(r.get_id(), "run1")
        self.assertEqual(r.get_openswath_filename(), "file1.csv")
        self.assertEqual(r.get_aligned_filename(), "file1.tsv")

    def test_get_peptide(self):
        r = Run([], {}, "run1", "file1.txt", filename="file1.csv", aligned_filename="file1.tsv")
        r.all_precursor_groups_ = dict( [ (str(i), MockPrecursorGroup(i)) for i in range(5) ]  )
        self.assertEqual( r.getPrecursorGroup("2").id_, 2) 
        self.assertIsNone( r.getPrecursorGroup("9_dummy"))

        ids = sorted([p.id_ for p in r])
        self.assertEqual( ids, list(range(5)))

    def test_get_best_peaks(self):
        r = Run([], {}, "run1", "file1.txt", filename="file1.csv", aligned_filename="file1.tsv")
        r.all_precursor_groups_ = dict( [ (str(i), MockPrecursorGroup(i)) for i in range(5) ]  )
        self.assertEqual( r.get_best_peaks(), ["42" for i in range(5)] )

if __name__ == '__main__':
    unittest.main()

