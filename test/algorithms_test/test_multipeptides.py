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

from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide

class MockPeakGroup():

    def __init__(self, fdr=1):
        self.fdr_score = fdr

    def get_fdr_score(self):
        return self.fdr_score

    def get_decoy(self):
        return False

class MockPeptide():

    def __init__(self, peakgroups, sequence):
        self.peakgroups = peakgroups
        self.id = -1
        self.sequence = sequence

    def get_all_peakgroups(self):
        return self.peakgroups

    def add_peakgroup(self, p):
        return self.peakgroups.append(p)

    def get_best_peakgroup(self):
        return self.peakgroups[0]

    def get_selected_peakgroup(self):
        return self.peakgroups[0]

    def get_decoy(self):
        return False

    def get_id(self):
        return "144"

class MockPrecursorGroup():

    def __init__(self, peptides, id_="dummy"):
        self.peptides = peptides
        self.id_ = id_

    def getPeptideGroupLabel(self):
        return self.id_

    def __lt__(self, other):
        return self.id_ > other.id_

    def __iter__(self):
        for p in self.peptides:
            yield p

    def get_decoy(self):
        return False


def help_insert(m):
    peakgroup = MockPeakGroup(0.2)
    mockPeptide2 = MockPeptide([ peakgroup ], "PEPTIDE_seq2")
    mockPrecursorGroup2 = MockPrecursorGroup([ mockPeptide2 ])

    peakgroup = MockPeakGroup(0.1)
    mockPeptide3 = MockPeptide([peakgroup], "PEPTIDE_seq2")
    mockPrecursorGroup3 = MockPrecursorGroup([ mockPeptide3 ])

    m.insert("42_0", mockPrecursorGroup2)
    m.insert("43_0", mockPrecursorGroup3)

class TestMultiPeptide(unittest.TestCase):

    def setUp(self):
        peakgroups = [MockPeakGroup() for i in range(2)]
        self.mockPeptide = MockPeptide(peakgroups, "PEPTIDESEQ")
        self.mockPrecursorGroup = MockPrecursorGroup([self.mockPeptide], "gr1")

    def testNrRuns(self):
        m = Multipeptide()
        m.set_nr_runs(42)
        self.assertEqual(m.get_nr_runs(), 42)

    def test_str(self):
        m = Multipeptide()
        myS = str(m)
        self.assertTrue(True)

    def test_getPrecursorGroups(self):
        m = Multipeptide()
        self.assertEqual(len(m.getPrecursorGroups()), 0)

    def test_insert(self):
        m = Multipeptide()
        self.assertEqual(len(m.getPrecursorGroups()), 0)

        m.insert("42_0", self.mockPrecursorGroup)
        self.assertEqual(m.getPrecursorGroup("42_0"), self.mockPrecursorGroup)
        self.assertEqual(len(list(m.getPrecursorGroup("42_0"))), 1)
        self.assertEqual(len( list(m.getPrecursorGroup("42_0"))[0].get_all_peakgroups() ), 2)
        self.assertEqual(len(m.getPrecursorGroups()), 1)
        self.assertTrue(m.hasPrecursorGroup("42_0"))

        # try to add more peakgroups to an already existing run
        peakgroups = [MockPeakGroup() for i in range(3)]
        mockPeptide2 = MockPeptide(peakgroups, "pepseq2")
        mockPrecursorGroup2 = MockPrecursorGroup([mockPeptide2], "gr2")

        myS = str(m)
        self.assertTrue(True)

    def test_insert_None(self):
        m = Multipeptide()
        self.assertEqual(len(m.getPrecursorGroups()), 0)

        m.insert("42_0", None)
        self.assertTrue(m.has_null_peptides())
        self.assertFalse(m.hasPrecursorGroup("42_0"))

    def test_getId(self):
        m = Multipeptide()
        self.assertIsNone(m.get_id())
        m.insert("42_0", self.mockPrecursorGroup)
        self.assertEqual(m.get_id(), "144")

    def test_more_than_fraction_selected(self):
        m = Multipeptide()
        self.assertIsNone(m.get_id())
        m.insert("42_0", self.mockPrecursorGroup)
        m.set_nr_runs(1)
        self.assertTrue(m.more_than_fraction_selected(0.1))
        self.assertTrue(m.more_than_fraction_selected(0.6))
        m.set_nr_runs(2)
        self.assertTrue(m.more_than_fraction_selected(0.1))
        self.assertFalse(m.more_than_fraction_selected(0.6))

    def test_get_decoy(self):
        m = Multipeptide()
        self.assertFalse(m.get_decoy())
        m.insert("42_0", self.mockPrecursorGroup)
        self.assertFalse(m.get_decoy())

    def test_all_above_cutoff(self):
        m = Multipeptide()
        m.set_nr_runs(2)
        self.assertFalse(m.all_above_cutoff(0.4))
        help_insert(m)

        self.assertTrue(m.all_above_cutoff(0.4))
        self.assertFalse(m.all_above_cutoff(0.15))

    def test_find_best_peptide_pg(self):

        m = Multipeptide()
        m.set_nr_runs(2)
        self.assertIsNone(m.find_best_peptide_pg())
        help_insert(m)

        self.assertAlmostEqual(m.find_best_peptide_pg().get_fdr_score(), 0.1)

    def test_all_selected(self):
        m = Multipeptide()
        m.set_nr_runs(2)
        self.assertIsNone(m.find_best_peptide_pg())

        help_insert(m)
        self.assertTrue(m.all_selected())

        m.set_nr_runs(3)
        self.assertFalse(m.all_selected())

if __name__ == '__main__':
    unittest.main()
