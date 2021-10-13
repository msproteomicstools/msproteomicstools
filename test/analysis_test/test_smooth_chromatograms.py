#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta test_smooth_chromatograms.py>
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
from analysis.chromatogram_utils.smooth_chromatograms import chromSmoother
import numpy as np

def Almost_equal_XICs(self, xic1, xic2, decimal = 6):
    for i in range(len(xic2)):
        for j in range(len(xic2[i][0])):
            self.assertAlmostEqual(xic2[i][0][j], xic1[i][0][j],  places = decimal)
            self.assertAlmostEqual(xic2[i][1][j], xic1[i][1][j],  places = decimal)

class TestSmoothXICs(unittest.TestCase):
    def setUp(self):
        x = np.arange(3003.4, 3048, 3.4)
        y = np.array([0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
                 4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923])
        self.XICs = [(x,y), (x,y), (x,y)]
    
    def test_smoothXICs(self):
        sm = chromSmoother(smoother = "sgolay", kernelLen = 9, polyOrd = 5)
        XICs_sm = sm.smoothXICs(self.XICs)
        XIC_true = (self.XICs[0][0], np.array([0.20636662, 0.88411087, 2.19019973, 3.76695006,
                          5.12687085, 5.77230554, 5.56200672, 4.5968725 , 3.2886408 , 1.97239146,
                          0.93076564, 0.34700936, 0.19229358, 0.14383756]))
        XICs_true = [XIC_true, XIC_true, XIC_true]
        Almost_equal_XICs(self, XICs_sm, XICs_true, decimal = 7)

if __name__ == '__main__':
    unittest.main()