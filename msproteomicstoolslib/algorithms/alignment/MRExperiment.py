#!/usr/bin/env python
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

from sys import stdout

import csv, os
import numpy
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.format.SWATHScoringReader import Run

class MRExperiment(object):
    """
    An MR (multirun) Experiment is a container for multiple experimental runs -
    some of which may contain the same precursors.

    # Read the files
    >>> fdr_cutoff = 0.01 # 1% FDR
    >>> reader = SWATHScoringReader.newReader(infiles, "openswath")
    >>> this_exp = Experiment()
    >>> this_exp.set_runs( reader.parse_files(options.realign_runs) )
    >>> multipeptides = this_exp.get_all_multipeptides(fdr_cutoff)
    """

    def __init__(self):
        self.runs = []

    def set_runs(self, runs):
        """Initialize with a set of runs.

        Args:
            runs(list(SWATHScoringReader.Run))
        """
        self.runs = runs

    def get_all_multipeptides(self, fdr_cutoff, verbose=False, verbosity=0):
        """Match all precursors in different runs to each other.

        Find all precursors that are above the fdr cutoff in each run and build
        a union of those precursors. Then search for each of those precursors
        in all the other runs and build a multipeptide / multiprecursor.
        """
        union_transition_groups = []
        union_proteins = []
        union_target_transition_groups = []
        for i,r in enumerate(self.runs):
            if verbose or verbosity >= 10: 
                stdout.write("\rParsing run %s out of %s" % (i+1, len(self.runs) ))
                stdout.flush()
            union_target_transition_groups.append( [peak.peptide.get_id() for peak in r.get_best_peaks_with_cutoff(fdr_cutoff) if not peak.peptide.get_decoy()] )
            union_transition_groups.append( [peak.peptide.get_id() for peak in r.get_best_peaks_with_cutoff(fdr_cutoff)] )
            union_proteins.append( list(set([peak.peptide.protein_name for peak in r.get_best_peaks_with_cutoff(fdr_cutoff) if not peak.peptide.get_decoy()])) )
        if verbose or verbosity >= 10: stdout.write("\r\r\n") # clean up

        union_target_transition_groups_set = set(union_target_transition_groups[0])
        self.union_transition_groups_set = set(union_transition_groups[0])
        self.union_proteins_set = set(union_proteins[0])
        for groups in union_transition_groups:
          self.union_transition_groups_set = self.union_transition_groups_set.union( groups )
        for groups in union_target_transition_groups:
          union_target_transition_groups_set = union_target_transition_groups_set.union( groups )
        for proteins in union_proteins:
          self.union_proteins_set = self.union_proteins_set.union( proteins )

        all_prec = sum([len(s) for s in union_transition_groups])
        target_prec = sum([len(s) for s in union_target_transition_groups])

        if verbose or verbosity >= 1:
            print "==================================="
            print "Finished parsing, number of precursors and peptides per run"
            print "All precursors", [len(s) for s in union_transition_groups], "(union of all runs %s)" % len(self.union_transition_groups_set)
            print "All target precursors", [len(s) for s in union_target_transition_groups], "(union of all runs %s)" % len(union_target_transition_groups_set)
            print "All target proteins", [len(s) for s in union_proteins], "(union of all runs %s)" % len(self.union_proteins_set)
            if all_prec > 0:
                print "Decoy percentage on precursor level %0.4f%%" % ( (all_prec - target_prec) * 100.0 / all_prec )

        self.initial_fdr_cutoff = fdr_cutoff
        if all_prec > 0 and all_prec - target_prec != 0:
            self.estimated_decoy_pcnt =  (all_prec - target_prec) * 100.0 / all_prec 
        else:
            self.estimated_decoy_pcnt = None

        multipeptides = []
        for peptide_id in self.union_transition_groups_set:
            m = Multipeptide()
            for r in self.runs:
                peptide = r.get_peptide(peptide_id)
                m.insert(r.get_id(), peptide)
            m.set_nr_runs(len(self.runs))
            multipeptides.append(m)
        return multipeptides


