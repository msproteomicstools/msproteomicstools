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
from nose.plugins.attrib import attr

import msproteomicstoolslib.math.Smoothing as smoothing

class TestUnitSmoothing(unittest.TestCase):

    def setUp(self):
        self.data1 = [5,7,8,9,10,15,7.1,6]
        self.data2 = [4,7,9,11,11,14,7.1,6.5]

    @attr('slow')
    def test_smooth_spline_scikit(self):
        """Test the smoothing spline using scikit"""
        sm = smoothing.SmoothingPy()
        r = sm._smooth_spline_scikit(self.data1, self.data2)

        self.assertEqual(len(r), 8)
        self.assertAlmostEqual(r[0], 4.41118020)
        self.assertAlmostEqual(r[2], 8.7900826361)
        self.assertAlmostEqual(r[5], 14.1276411901)
        self.assertAlmostEqual(r[7], 5.90758396)

        r = sm._smooth_spline_scikit(self.data1, self.data2, [5,7,10.0])
        self.assertEqual(len(r), 3)
        self.assertAlmostEqual(r[0], 4.6129201633)
        self.assertAlmostEqual(r[1], 7.73837621136)
        self.assertAlmostEqual(r[2], 10.3726686328)

        r = sm._smooth_spline_scikit(self.data1, self.data2, [10.0,5.0,7])
        self.assertEqual(len(r), 3)
        self.assertAlmostEqual(r[0], 10.3726686328)
        self.assertAlmostEqual(r[1], 4.6129201633)
        self.assertAlmostEqual(r[2], 7.73837621136)

        # Since the (xhat,r) is optimized, different lamda values are estimated
        # and a different estimation is used here...
        r = sm._smooth_spline_scikit(self.data1, self.data2, [10.0,5.0,7,10.00001])
        self.assertEqual(len(r), 4)
        self.assertAlmostEqual(r[0], 11.60276344181, 3)
        self.assertAlmostEqual(r[1], 4.327929410625, 3)
        self.assertAlmostEqual(r[2], 7.360352947814, 3)
        self.assertAlmostEqual(r[3], 11.60276823628, 3)

        r = sm._smooth_spline_scikit(self.data1, self.data2, [10.0,5.0,7,10.0], True)
        self.assertEqual(len(r), 4)
        expected = [10.372668632871067, 4.6129201633602159, 7.738376211369804, 10.372668632871067]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res, exp, 4)

    @attr('slow')
    def test_smooth_spline_scikit_wrap(self):
        """Test the smoothing spline using scikit and the wrapper"""
        sm = smoothing.SmoothingPy()
        import numpy
        sm.initialize(self.data1, self.data2, xmin=numpy.min(numpy.array(self.data1)), xmax=numpy.max(numpy.array(self.data1)))
        r_pred = sm.predict(self.data1)
        r = sm._smooth_scikit_legacy(self.data1, self.data2, self.data1)
        self.assertEqual(len(r), 8)
        expected = [4.340432925607892, 7.412884078435387, 8.865581929053054,
                    10.157652480992748, 11.196560644581082, 14.158931578788266,
                    7.561207213410677, 5.90694677222806]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res,exp)
        for res, exp in zip(r_pred,expected):
            self.assertAlmostEqual(res,exp)

        r = sm._smooth_scikit_legacy(self.data1, self.data2, [5,7,10.0])
        self.assertEqual(len(r), 3)
        expected = [4.307319862409416, 7.423679061650855, 11.15630836580813]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res,exp)

        r = sm._smooth_scikit_legacy(self.data1, self.data2, [10.0,5.0,7])
        self.assertEqual(len(r), 3)
        expected = [11.15630836580813, 4.307319862409416, 7.423679061650855]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res,exp)

        # Here, we expect the exact same results
        r = sm._smooth_scikit_legacy(self.data1, self.data2, [10.0,5.0,7,10.0])
        self.assertEqual(len(r), 4)
        expected = [11.15630836580813, 4.307319862409416, 7.423679061650855, 11.15630836580813]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res,exp)

        # However, if we also evaluate at a point outside the previous window, we expect the results to change slightly
        r = sm._smooth_scikit_legacy(self.data1, self.data2, [10.0,5.0,7,10.0, 15.0])
        self.assertEqual(len(r), 5)
        expected = [11.196560644581082, 4.340432925607892, 7.412884078435387, 11.196560644581082, 14.158931578788266]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res,exp)

        # If we chose a point that is very far away, we still expect a "reasonable" result (e.g. -100 becores -96)
        r = sm._smooth_scikit_legacy(self.data1, self.data2, [10.0,5.0,7,-100.0, 15.0])
        self.assertEqual(len(r), 5)
        expected = [10.265638884711272, 5.411029040286351, 7.352910860216361, -96.53810165019972, 15.119857967104132]
        for res, exp in zip(r,expected):
            self.assertAlmostEqual(res,exp, 1)

    def test_smooth_spline_scipy_uni(self):
        """Test the univariate spline using spline (no crossvalidation)"""
        sm = smoothing.UnivarSplineNoCV()
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data1)
        self.assertEqual(len(r), 8)

        expected =  [4.1457135130652265, 7.442333705513172, 8.9177273726462474, 10.248891199849604, 11.414111021721407, 13.991054262042756, 7.5957445613093642, 5.8444243638522186]
        self.assertEqual(len(r), 8)
        for a,b in zip(expected,r):
          self.assertAlmostEqual(a, b)

        r = sm.predict(self.data2)
        expected = [2.34266247,   7.2926131 ,  10.48943975,  11.85840597,
                11.85840597,  13.48225519,   7.44184246,   6.61579704]

        self.assertEqual(len(r), 8)
        for a,b in zip(expected,r):
          self.assertTrue( abs(1-a/b) < 0.05)

    def test_smooth_spline_scipy_cv(self):
        """Test the univariate spline using spline (with crossvalidation)"""
        sm = smoothing.UnivarSplineCV()
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data1)
        self.assertEqual(len(r), 8)

        # This is slightly better than the NoCV variation
        expected =  [4.1626385843797094, 7.3804099239612029, 8.9667396489152544, 10.384851777100122, 11.316311505414465, 13.994282700490476, 7.5367306411050095, 5.8580352186337672]
        for a,b in zip(expected,r):
            self.assertTrue( abs(1.0-a/b) < 0.1)

    def test_smooth_nn(self):
        """Test the univariate spline using local kernel smoothing"""

        # Test with regular parameters
        sm = smoothing.WeightedNearestNeighbour(2, 5, 0.5, False)
        sm.initialize(self.data1, self.data2)
        r = sm.predict( [ 5 ] )
        self.assertAlmostEqual(r[0], 4.85378590078)

        r = sm.predict( [ 15 ] )
        self.assertAlmostEqual(r[0], 14.472485768500951)

        # Test with exponent
        sm = smoothing.WeightedNearestNeighbour(2, 5, 0.5, False, exponent=2.0)
        sm.initialize(self.data1, self.data2)
        r = sm.predict( [ 5 ] )
        self.assertAlmostEqual(r[0], 4.4223582231809182)

        r = sm.predict( [ 15 ] )
        self.assertAlmostEqual(r[0], 14.04993649085635)

        sm = smoothing.WeightedNearestNeighbour(3, 5, 0.5, False)
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data1)
        self.assertEqual(len(r), 8)

        # This is slightly better than the NoCV variation
        expected = [4.85378590078329, 7.3181818181818183, 8.6853448275862046, 10.054730258014073, 11.044451654769629, 14.497816294331514, 7.4375518352136076, 6.2364096080910238]

        for a,b in zip(expected,r):
            self.assertAlmostEqual(a, b)

        # Try with smaller mindiff => more weight on close neighbors
        sm = smoothing.WeightedNearestNeighbour(3, 5, 0.1, False)
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data1)
        self.assertEqual(len(r), 8)

        # This is slightly better than the NoCV variation
        expected = [4.3099526066350711, 7.0999999999999996, 8.8596153846153847, 10.610377054463424, 11.015838150289017, 14.1233812970026, 7.2064265084789056, 6.3871142393069835]
        for a,b in zip(expected,r):
            self.assertAlmostEqual(a, b)

    def test_smooth_lowess(self):
        """Test the lowess smoothing"""
        sm = smoothing.LowessSmoothingBiostats()
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data1)
        self.assertEqual(len(r), 8)

        # This is slightly better than the NoCV variation
        expected = [4.2123769729879061, 7.3305230831706876, 8.8162015867770727, 10.144542883530072, 11.507036080352814, 14.061195393451431, 7.4880128821482463, 5.7476045207327786]
        for a,b in zip(expected,r):
            self.assertAlmostEqual(a, b)

    @attr('rpy2')
    def test_smooth_spline_r(self):
        """Test the smoothing spline using R"""
        sm = smoothing.SmoothingR()
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data2)
        expected = [  2.34266247,   7.2926131 ,  10.48943975,  11.85840597,
                11.85840597,  13.48225519,   7.44184246,   6.61579704]

        self.assertEqual(len(r), 8)
        for a,b in zip(expected,r):
          self.assertAlmostEqual(a, b)

    def test_smooth_spline_r_extern(self):
        """Test the smoothing spline using external R"""
        sm = smoothing.SmoothingRExtern()
        sm.initialize(self.data1, self.data2)
        r = sm.predict(self.data2)
        expected = [  2.34266247,   7.2926131 ,  10.48943975,  11.85840597,
                11.85840597,  13.48225519,   7.44184246,   6.61579704]

        self.assertEqual(len(r), 8)
        for a,b in zip(expected,r):
          self.assertAlmostEqual(a, b)

    def test_duplication(self):
        """
        Test de-duplication of array data
        """
        arr = [0, 0, 5, 6, 6, 7, 8, 8]
        sm = smoothing.SmoothingPy()
        de_dupl,duplications = sm.de_duplicate_array(arr)
        re_dupl = sm.re_duplicate_array(de_dupl, duplications)
        # the input and output need to be identical!
        self.assertEqual(re_dupl, arr)

    def test_gettingOperator(self):
        """
        Test getting the correct smoothing operator
        """
        op = smoothing.get_smooting_operator(use_linear=True)
        self.assertTrue(isinstance(op, smoothing.SmoothingLinear))

        op = smoothing.get_smooting_operator(use_scikit=True)
        self.assertTrue(isinstance(op, smoothing.SmoothingPy))

        op = smoothing.get_smooting_operator(use_external_r=True, tmpdir="tmp")
        self.assertTrue(isinstance(op, smoothing.SmoothingRExtern))

    @attr('rpy2')
    def test_gettingOperator_rpy2(self):
        """
        Test getting the correct smoothing operator
        """
        op = smoothing.get_smooting_operator()
        self.assertTrue(isinstance(op, smoothing.SmoothingR))

        op = smoothing.getSmoothingObj("splineR")
        self.assertTrue(isinstance(op, smoothing.SmoothingR))

    def test_gettingOperator_obj(self):
        """
        Test getting the correct smoothing operator (new interface)
        """

        op = smoothing.getSmoothingObj("diRT")
        self.assertTrue(isinstance(op, smoothing.SmoothingNull))

        op = smoothing.getSmoothingObj("None")
        self.assertTrue(isinstance(op, smoothing.SmoothingNull))

        op = smoothing.getSmoothingObj("linear")
        self.assertTrue(isinstance(op, smoothing.SmoothingLinear))

        # op = smoothing.getSmoothingObj("splineR")
        # self.assertTrue(isinstance(op, smoothing.SmoothingR))

if __name__ == '__main__':
    unittest.main()
