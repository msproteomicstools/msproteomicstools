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

class TestDAMReader(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "msproteomicstoolslib/format/")

    def exact_diff(self, name1, name2):
        f1 = open(name1, "r")
        f2 = open(name2, "r")
        for l1,l2 in zip(f1,f2):
            for field1,field2 in zip(l1.split(),l2.split()):
                try:
                    self.assertAlmostEqual(float(field1),float(field2) )
                except ValueError:
                    self.assertEqual(field1,field2)

    def test_damReader(self):
        script = os.path.join(os.path.join(self.scriptdir, ""), "methodDamReader.py")
        filename = os.path.join(self.datadir, "ABSciex_testInput.dam")
        expected_outcome = os.path.join(self.datadir, "methodDamReader_1_output.csv")
        tmpfilename = "methodfile.out.tmp"

        args = "--doAssert --in %s --out %s" % (filename, tmpfilename)
        cmd = "python %s %s" % (script, args)
        sub.check_call(cmd,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
        
        self.exact_diff(tmpfilename, expected_outcome)

        os.remove(tmpfilename)

    def test_damReader2(self):
        #tests for the correct reading after a large "FFFFF" insertion into the method parameters
        script = os.path.join(os.path.join(self.scriptdir, ""), "methodDamReader.py")
        filename = os.path.join(self.datadir, "ABSciex_testInput2.dam")
        expected_outcome = os.path.join(self.datadir, "methodDamReader_2_output.csv")
        tmpfilename = "methodfile.out.tmp"

        args = "--doAssert --in %s --out %s" % (filename, tmpfilename)
        cmd = "python %s %s" % (script, args)
        sub.check_call(cmd,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
        
        self.exact_diff(tmpfilename, expected_outcome)

        os.remove(tmpfilename)

    def test_damReader3(self):
        #tests for the correct translation of a non-ASCII character ('ð') in the name/ID of the method paramters
        script = os.path.join(os.path.join(self.scriptdir, ""), "methodDamReader.py")
        filename = os.path.join(self.datadir, "ABSciex_testInput3.dam")
        expected_outcome = os.path.join(self.datadir, "methodDamReader_3_output.csv")
        tmpfilename = "methodfile.out.tmp"

        args = "--doAssert --in %s --out %s" % (filename, tmpfilename)
        cmd = "python %s %s" % (script, args)
        sub.check_call(cmd,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
        
        self.exact_diff(tmpfilename, expected_outcome)

        os.remove(tmpfilename)

if __name__ == '__main__':
    unittest.main()
