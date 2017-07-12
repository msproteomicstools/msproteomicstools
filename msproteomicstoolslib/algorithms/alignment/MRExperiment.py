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

from __future__ import print_function
from sys import stdout

import csv, os
import numpy
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.format.SWATHScoringReader import Run

class MRExperiment(object):
    """
    An MR (multirun) Experiment is a container for multiple experimental runs.

    In some of the runs the same peptidde precursors may be identified and the
    job of this object is to keep track of these experiments and the identified
    precursors across multiple runs.

    Example usage:

    >>> # Read the files
    >>> fdr_cutoff = 0.01 # 1% FDR
    >>> reader = SWATHScoringReader.newReader(infiles, "openswath")
    >>> this_exp = Experiment()
    >>> this_exp.set_runs( reader.parse_files(options.realign_runs) )
    >>> multipeptides = this_exp.get_all_multipeptides(fdr_cutoff)
    """

    def __init__(self):
        self.runs = []
        self.nr_ambiguous = -1
        self.nr_multiple_align = -1

    def set_runs(self, runs):
        """Initialize with a set of runs.

        Args:
            runs(list of :class:`.Run`) : A set of runs
        """
        self.runs = runs

    def get_all_multipeptides(self, fdr_cutoff, verbose=False, verbosity=0):
        """Match all precursors in different runs to each other.

        Find all precursors that are above the fdr cutoff in each run and build
        a union of those precursors. Then search for each of those precursors
        in all the other runs and build a multipeptide / multiprecursor.

        Parameters
        ----------
        fdr_cutoff : float
            A cutoff in fdr (between 0 and 1) to use for the alignment. Each
            generated Multipeptide needs to have at least one member who is below
            the cutoff.
        verbose : bool
            Whether to be verbose or not
        verbosity : int
            How verbose to be
        """

        # Identify across all runs which transition groups are above the cutoff
        union_transition_groups = []
        union_target_transition_groups = []
        union_proteins = []

        self.union_transition_groups_set = set([])
        self.union_proteins_set = set([])
        self.union_target_transition_groups_set = set()
        for i,r in enumerate(self.runs):
            gr = []
            gr_target = []
            gr_protein = []
            for precursor_group in r:
                for peptide_precursor in precursor_group:
                    if (peptide_precursor.get_best_peakgroup().get_fdr_score() < fdr_cutoff):
                        gr.append( precursor_group.getPeptideGroupLabel() )
                        if not precursor_group.get_decoy():
                            gr_target.append(precursor_group.getPeptideGroupLabel())
                            gr_protein.append(peptide_precursor.getProteinName())
            union_transition_groups.append(gr)
            union_target_transition_groups.append(gr_target)
            union_proteins.append(list(set(gr_protein)))

            self.union_target_transition_groups_set = self.union_target_transition_groups_set.union(gr_target)
            self.union_transition_groups_set = self.union_transition_groups_set.union(gr)
            self.union_proteins_set = self.union_proteins_set.union(gr_protein)

        if verbose or verbosity >= 10: 
            stdout.write("\r\r\n") # clean up

        all_prec = sum([len(s) for s in union_transition_groups])
        target_prec = sum([len(s) for s in union_target_transition_groups])

        if verbose or verbosity >= 1:
            print("===================================")
            print("Finished parsing, number of precursors and peptides per run")
            print("All precursors", [len(s) for s in union_transition_groups], "(union of all runs %s)" % len(self.union_transition_groups_set))
            print("All target precursors", [len(s) for s in union_target_transition_groups], "(union of all runs %s)" % len(self.union_target_transition_groups_set))
            print("All target proteins", [len(s) for s in union_proteins], "(union of all runs %s)" % len(self.union_proteins_set))
            if all_prec > 0:
                print("Decoy percentage on precursor level %0.4f%%" % ( (all_prec - target_prec) * 100.0 / all_prec ))

        self.initial_fdr_cutoff = fdr_cutoff
        if all_prec > 0 and all_prec - target_prec != 0:
            self.estimated_decoy_pcnt =  (all_prec - target_prec) * 100.0 / all_prec 
        else:
            self.estimated_decoy_pcnt = None

        multipeptides = []
        for peptide_id in self.union_transition_groups_set:
            m = Multipeptide()
            for r in self.runs:
                precursor_group = r.getPrecursorGroup(peptide_id)
                m.insert(r.get_id(), precursor_group)
            m.set_nr_runs(len(self.runs))
            multipeptides.append(m)

        # Return sorted multipeptides for consistency across all Python versions
        return(sorted(multipeptides, key=lambda x: str(x)))

