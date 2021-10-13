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
import os
from analysis.alignment.reference_run_selection import referenceForPrecursor
from msproteomicstoolslib.data_structures.Run import Run

class TestReferenceRunSelection(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader
        cls.dirname = os.path.dirname(os.path.abspath(__file__))
        cls.topdir = os.path.join(os.path.join(cls.dirname, ".."), "..")
        cls.datadir = os.path.join(os.path.join(cls.topdir, "test"), "data")
        cls.datadir_DIAlign = os.path.join(cls.datadir, "DIAlign")

        filename = os.path.join(cls.datadir_DIAlign, "merged.osw")
        r = SWATHScoringReader.newReader([filename], "openswath", "minimal")
        runs = r.parse_files(read_exp_RT = False)
        from analysis.alignment.feature_alignment import Experiment        
        this_exp = Experiment()
        this_exp.set_runs(runs)
        cls.best_run = this_exp.determine_best_run(alignment_fdr_threshold = 0.05)
        
        cls.mp = this_exp.get_all_multipeptides(0.05, verbose=False)

    def test_get_reference_for_precursors(self):

        prec_ref = referenceForPrecursor("precursor_specific", self.best_run, 0.05)
        reference_run = prec_ref.get_reference_for_precursors(self.mp)
        self.assertEqual(prec_ref.best_run.get_id(), self.best_run.get_id())
        self.assertEqual(len(reference_run), 202)
        
        reference_run = referenceForPrecursor("precursor_specific", None, 0.05).get_reference_for_precursors(self.mp)
        self.assertEqual(len(reference_run), 202)
        self.assertIsInstance(reference_run[32], Run)
        self.assertEqual(reference_run[32].get_id(), 2234664662238281994)

        reference_run = referenceForPrecursor("multipeptide_specific", None, 0.01).get_reference_for_precursors(self.mp)
        self.assertEqual(len(reference_run), 202)
        self.assertIsNone(reference_run[32])
        # print(reference_run)
        # self.assertEqual(reference_run[32].get_id(), 2234664662238281994)

if __name__ == '__main__':
    unittest.main()