#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta reference_run_selection.py>
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
from msproteomicstoolslib.util import getBaseName

class TestBaseName(unittest.TestCase):
    def test_getBaseName(self):
        filename = 'data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz'
        self.assertEqual(getBaseName(filename), 'hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt')

        filename = 'data\\raw\\hr_K120808_Strep10%P.sqMass'
        self.assertEqual(getBaseName(filename), 'hr_K120808_Strep10%P')

if __name__ == '__main__':
    unittest.main()