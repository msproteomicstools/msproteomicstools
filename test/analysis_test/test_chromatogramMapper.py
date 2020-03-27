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
from analysis.chromatogram_utils.chromatogramMapper import mz_pointer, mzml_accessors
from msproteomicstoolslib.data_structures.Run import Run

def Almost_equal_XICs(self, xic1, xic2, decimal = 6):
    for i in range(len(xic2)):
        for j in range(len(xic2[i][0])):
            self.assertAlmostEqual(xic2[i][0][j], xic1[i][0][j],  places = decimal)
            self.assertAlmostEqual(xic2[i][0][j], xic1[i][0][j],  places = decimal)

def get_XICs(self):
    import pickle
    with open(os.path.join(self.datadir_DIAlign,'XIC_group_run2.data'), 'rb') as filehandle:
        XIC_group_run2 = pickle.load(filehandle)
    return XIC_group_run2


class TestMzPointer(unittest.TestCase):
    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.datadir_DIAlign = os.path.join(self.datadir, "DIAlign") # Instance attribute
    
    def test_init(self):
        filename = os.path.join(self.datadir_DIAlign, 'hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML')
        mz = mz_pointer(filename)
        self.assertEqual(mz.filename, filename)
        self.assertEqual(mz.maxIndex, 71)
    
    def test_extractXIC_group_(self):
        filename = os.path.join(self.datadir_DIAlign, 'hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML')
        mz = mz_pointer(filename)
        # Precursor = 4618
        # TransitionID = [27706, 27707, 27708, 27709, 27710, 27711]
        chromIndices = [36, 37, 38, 39, 40, 41]
        XICs = mz.extractXIC_group_(chromIndices)
        """
        import pickle
        XIC_group_run2 = mz.extractXIC_group_(chromIndices)
        with open(os.path.join(self.datadir_DIAlign,'XIC_group_run2.data'), 'wb') as filehandle:
            pickle.dump(XIC_group_run2, filehandle)
        """
        XIC_group_run2 = get_XICs(self)
        Almost_equal_XICs(self, XICs, XIC_group_run2, decimal = 6)

class TestMzmlAccessor(unittest.TestCase):
    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.topdir = os.path.join(os.path.join(self.dirname, ".."), "..")
        self.datadir = os.path.join(os.path.join(self.topdir, "test"), "data")
        self.datadir_DIAlign = os.path.join(self.datadir, "DIAlign") # Instance attribute
        filename = os.path.join(self.datadir_DIAlign, 'merged.osw')
        self.chromFile0 = os.path.join(self.datadir_DIAlign, 'hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML')
        self.chromFile2 = os.path.join(self.datadir_DIAlign, 'hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML')
        run0 = Run([], {}, 125704171604355508, filename, 'data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz',
         'data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz', useCython=False)
        run2 = Run([], {}, 2234664662238281994, filename, 'data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz',
         'data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz', useCython=False)
        self.runs = [run0, run2]
        self.MStoFeature = {'data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz': (self.chromFile0, run0),
        'data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz':(self.chromFile2, run2)}

    def test_init(self):
        accsr = mzml_accessors(self.runs, self.MStoFeature)
        self.assertEqual(accsr.pointers[125704171604355508].filename, self.chromFile0)
        self.assertEqual(accsr.pointers[2234664662238281994].filename, self.chromFile2)
        self.assertEqual(accsr.runs, self.runs)
        self.assertEqual(accsr.run_chromIndex_map, {})

    def test_extractXIC_group(self):
        accsr = mzml_accessors(self.runs, self.MStoFeature)
        accsr.set_precursor_to_chromID({4618:[27706, 27707, 27708, 27709, 27710, 27711]})
        XICs = accsr.extractXIC_group(self.runs[1], 4618)
        XIC_group_run2 = get_XICs(self)
        Almost_equal_XICs(self, XICs, XIC_group_run2, decimal = 6)

    def test_set_precursor_to_chromID(self):
        accsr = mzml_accessors(self.runs, self.MStoFeature)
        accsr.set_precursor_to_chromID({4618:[27706, 27707, 27708, 27709, 27710, 27711]})
        accsr.run_chromIndex_map[2234664662238281994][4618]
        self.assertEqual(accsr.run_chromIndex_map[125704171604355508][4618], [36, 37, 38, 39, 40, 41])
        self.assertEqual(accsr.run_chromIndex_map[2234664662238281994][4618], [36, 37, 38, 39, 40, 41])

if __name__ == '__main__':
    unittest.main()