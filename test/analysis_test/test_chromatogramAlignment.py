#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta test_chromatogramAlignment.py>
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
import subprocess as sub
from analysis.alignment import chromatogram_alignment as chromAlign
from msproteomicstoolslib.data_structures.Run import Run
from msproteomicstoolslib.data_structures.Precursor import Precursor
import numpy as np

class TestChromatogramAlignment(unittest.TestCase):
    def setUp(self):
        self.trgr_id = 32
        self.peptide_group_label = 4732
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.scriptdir = os.path.join(self.topdir, "analysis")
        self.datadir_DIAlign = os.path.join(self.datadir, "DIAlign") # Instance attribute
            
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
                    try:
                        self.assertAlmostEqual(float(field1),float(field2) )
                    except ValueError:
                        self.assertEqual(field1,field2)
                except AssertionError as e:
                    print("Assertion failed", field1, field2, "in line", l1, l2)
                    print("See diff between files: diff", name1, name2)
                    raise e
    
    def test_get_aligned_time(self):
        t = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
        indices = [-1, 0, 1, 2, -1, 3, 4, -1, 5, 6, 7, -1]
        t_aligned = chromAlign.get_aligned_time(t, indices, skipValue = -1)
        t_true = np.array([np.nan, 1.1, 1.2, 1.3, 1.35, 1.4, 1.5, 1.55, 1.6, 1.7, 1.8, np.nan])
        np.testing.assert_array_equal(t_true, t_aligned)

    def test_updateRetentionTime(self):
        run = Run([], {}, 125704171604355508, 'merged.osw', 'file.mzML.gz', 'file.mzML.gz', useCython=False)
        p = Precursor(self.trgr_id, run)
        run.addPrecursor(p, self.peptide_group_label)
        run.getPrecursor(self.peptide_group_label, self.trgr_id).add_peakgroup_tpl((364283, 0.001, 1.47, 3000), self.trgr_id, cluster_id = -1)
        t_ref = np.array([np.nan, 21.1, 21.2, 21.3, 21.35, 21.4, 21.5, 21.55, 21.6, 21.7, 21.8, np.nan])
        t_eXp = np.array([np.nan, 1.1, 1.2, 1.3, 1.35, 1.4, 1.5, 1.55, 1.6, 1.7, 1.8, np.nan])
        chromAlign.updateRetentionTime(run, self.peptide_group_label, self.trgr_id, t_ref, t_eXp)
        self.assertEqual(run.getPrecursor(self.peptide_group_label, self.trgr_id).peakgroups_, [(364283, 0.001, 21.5, 3000, None)])

    def test_1_chromatogram_alignment(self):
        script = os.path.join(os.path.join(self.scriptdir, "alignment"), "chromatogram_alignment.py")
        expected_outcome = os.path.join(self.datadir_DIAlign, "chromatogram_alignment_1.csv")
        chrom1 = os.path.join(self.datadir_DIAlign, "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML")
        chrom2 = os.path.join(self.datadir_DIAlign, "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
        chrom3 = os.path.join(self.datadir_DIAlign, "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
        chroms = ' '.join([chrom1, chrom2, chrom3])
        featureFile = os.path.join(self.datadir_DIAlign, "merged.osw")
        tmpfilename = "chromatogram_alignment_1.tmp.csv"
        
        args = "--chromatograms %s --features %s --out_matrix %s" % (chroms, featureFile, tmpfilename)
        cmd = "python %s %s" % (script, args)
        
        sub.check_output(cmd, shell=True)
        self.exact_diff(tmpfilename, expected_outcome)
        
        os.remove(tmpfilename)


if __name__ == '__main__':
    unittest.main()
