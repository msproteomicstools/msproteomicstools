#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta test_chromatogramSmoothing.py>
Copyright (C) 2020 Shubham Gupta

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

--------------------------------------------------------------------------
$Maintainer: Shubham Gupta$
$Authors: Shubham Gupta$
--------------------------------------------------------------------------
"""
import unittest
from msproteomicstoolslib.math import ChromatogramSmoothing as smoother
import numpy as np

def Almost_equal_XIC(self, xic1, xic2, decimal = 6):
    for j in range(len(xic2[0])):
        self.assertAlmostEqual(xic2[0][j], xic1[0][j],  places = decimal)
        self.assertAlmostEqual(xic2[1][j], xic1[1][j],  places = decimal)


class TestXICsmoothig(unittest.TestCase):
    def setUp(self):
        x = np.arange(3003.4, 3048, 3.4)
        y = np.array([0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
                 4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923])
        self.XIC = (x,y)

    def test_LoessSmooth(self):
        sm = smoother.getXIC_SmoothingObj(smoother = "loess", kernelLen = 3.9, polyOrd = 1)
        sm.initialize(self.XIC[0], self.XIC[1])
        XIC_sm = sm.smooth(self.XIC[0], self.XIC[1])
        XIC_true = (self.XIC[0], np.array([0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                5.8288915, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                0.9520709, 0.3294218, 0.2009581, 0.1420923]))
        Almost_equal_XIC(self, XIC_sm, XIC_true, decimal = 7)

    def test_SgolaySmooth(self):
        sm = smoother.getXIC_SmoothingObj(smoother = "sgolay", kernelLen = 9, polyOrd = 5)
        sm.initialize(self.XIC[0], self.XIC[1])
        XIC_sm = sm.smooth(self.XIC[0], self.XIC[1])
        XIC_true = (self.XIC[0], np.array([0.20636662, 0.88411087, 2.19019973, 3.76695006,
                          5.12687085, 5.77230554, 5.56200672, 4.5968725 , 3.2886408 , 1.97239146,
                          0.93076564, 0.34700936, 0.19229358, 0.14383756]))
        Almost_equal_XIC(self, XIC_sm, XIC_true, decimal = 7)

    def test_GaussianSmooth(self):
        sm = smoother.getXIC_SmoothingObj(smoother = "gaussian", kernelLen = 4)
        sm.initialize(self.XIC[0], self.XIC[1])
        XIC_sm = sm.smooth(self.XIC[0], self.XIC[1])
        XIC_true = (self.XIC[0], np.array([0.29757935, 1.00811467, 2.24531529, 3.70565892, 5.01475393,
                        5.64511602, 5.40968113, 4.51346745, 3.29795671, 2.02014922,
                        1.02574116, 0.42391014, 0.21569467, 0.12735513]))
        Almost_equal_XIC(self, XIC_sm, XIC_true, decimal = 7)

    def test_UnknownMethod(self):
        self.assertRaises(Exception, smoother.getXIC_SmoothingObj, "Unknown")

if __name__ == '__main__':
    unittest.main()