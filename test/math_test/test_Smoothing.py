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

import msproteomicstoolslib.math.Smoothing as smoothing

class TestUnitSmoothing(unittest.TestCase):

    def setUp(self):
        self.data1 = [5,7,8,9,10,15,7.1,6]
        self.data2 = [4,7,9,11,11,14,7.1,6.5]

    def test_smooth_spline_scikit(self):
        """Test the smoothing spline using scikit"""
        sm = smoothing.Smoothing()
        r = sm.smooth_spline_scikit(self.data1, self.data2)

        self.assertEqual(len(r), 8)
        self.assertAlmostEqual(r[0], 4.41118020)
        self.assertAlmostEqual(r[2], 8.7900826361)
        self.assertAlmostEqual(r[5], 14.1276411901)
        self.assertAlmostEqual(r[7], 5.90758396)

        r = sm.smooth_spline(self.data1, self.data2, [5,7,10.0])
        self.assertEqual(len(r), 3)
        self.assertAlmostEqual(r[0], 4.6129201633)
        self.assertAlmostEqual(r[1], 7.73837621136)
        self.assertAlmostEqual(r[2], 10.3726686328)

        r = sm.smooth_spline(self.data1, self.data2, [10.0,5.0,7])
        self.assertEqual(len(r), 3)
        self.assertAlmostEqual(r[0], 10.3726686328)
        self.assertAlmostEqual(r[1], 4.6129201633)
        self.assertAlmostEqual(r[2], 7.73837621136)

        # Since the (xhat,r) is optimized, different lamda values are estimated
        # and a different estimation is used here...
        r = sm.smooth_spline(self.data1, self.data2, [10.0,5.0,7,10.00001])
        self.assertEqual(len(r), 4)
        self.assertAlmostEqual(r[0], 11.60276344181)
        self.assertAlmostEqual(r[1], 4.3279294106253)
        self.assertAlmostEqual(r[2], 7.3603529478143)
        self.assertAlmostEqual(r[3], 11.60276823628391)

        r = sm.smooth_spline(self.data1, self.data2, [10.0,5.0,7,10.0], True)
        self.assertEqual(len(r), 4)
        self.assertAlmostEqual(r[0], 11.60276344181)
        self.assertAlmostEqual(r[1], 4.3279294106253)
        self.assertAlmostEqual(r[2], 7.3603529478143)
        self.assertAlmostEqual(r[3], 11.60276823628391)

    def test_smooth_spline_r(self):
        """Test the smoothing spline using R"""
        sm = smoothing.Smoothing()
        r = sm.smooth_spline_r(self.data1, self.data2, self.data2)
        expected = [  2.34266247,   7.2926131 ,  10.48943975,  11.85840597,
                11.85840597,  13.48225519,   7.44184246,   6.61579704]

        self.assertEqual(len(r), 8)
        for a,b in zip(expected,r):
          self.assertAlmostEqual(a, b)

if __name__ == '__main__':
    unittest.main()
