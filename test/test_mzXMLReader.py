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

import msproteomicstoolslib.format.mzXMLreader as mzXMLReader

class TestmzXMLReader(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        joined = os.path.join(dirname, "data")
        # the test file is AFA1451A160MA2 of the dataset PAe001446_mzXML_201102101633.tar.gz from Peptide Atlas
        self.filename = os.path.join(joined, "testfile.small.mzXML")

    def test_readfile(self):
        reader = mzXMLReader.mzXMLReader(self.filename, True)

    def test_readscan(self):
        import msproteomicstoolslib.format.mzXMLreader as mzXMLReader
        reader = mzXMLReader.mzXMLReader(self.filename, True)
        scan = reader.read_scan(5)

    def test_readpeaks(self):
        import msproteomicstoolslib.format.mzXMLreader as mzXMLReader
        reader = mzXMLReader.mzXMLReader(self.filename, True)
        scan = reader.read_scan(5, True)
        self.assertEqual( len(scan.peaks), 615)
        self.assertAlmostEqual( scan.peaks[0].int, 271.413513184)
        self.assertAlmostEqual( scan.peaks[0].mz, 350.3074646)
        self.assertAlmostEqual( scan.max_peak().int, 10638.6044922)
        self.assertAlmostEqual( scan.max_peak().mz, 390.929901123)

    def test_parse_scans(self):
        reader = mzXMLReader.mzXMLReader(self.filename, False)
        scans = reader.parse_scans(ms2Only=False)

    def test_create_ms_rt_hashes(self):
        reader = mzXMLReader.mzXMLReader(self.filename, False)
        scans = reader.parse_scans()
        reader.create_ms_rt_hashes(scans)

if __name__ == '__main__':
    unittest.main()
