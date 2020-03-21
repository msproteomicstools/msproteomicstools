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

class referenceForPrecursor():
    """
    Calculates reference run for each precursor.
    Returns a dictionary with Precursor_id as key and Run as value.
    """
    def __init__(self, refType="best_run", run = None, alignment_fdr_threshold = 0.05):
        """
        refType must be either best_run, multipeptide_specific or precursor_specific. 
        """
        self.referenceType = refType
        self.alignment_fdr_threshold = alignment_fdr_threshold
        if(refType == "best_run"):
            if not run:
                raise Exception("run should be provided if refType is best_run.")
            self.best_run = run

    def get_reference_for_precursors(self, multipeptides, refType = "best_run"):
        if (refType == "best_run"):
            return self._get_reference_run(multipeptides)
        elif (refType == "precursor_specific"):
            return self._get_precursor_reference_run(multipeptides)
        elif (refType == "multipeptide_specific"):
            return self._get_multipeptide_reference_run(multipeptides)
        else:
            raise Exception("refType must be either best_run, multipeptide_specific or precursor_specific.")

    # Get a single reference run
    def _get_reference_run(self, multipeptides):
        """
        Returns a dectionary with precursor id as key and run as value.
        Single reference run for each precursor.
        """
        reference_run = {}
        run_id = self.best_run.get_id()
        precursor_ids = set()
        for i in range(len(multipeptides)):
            # Get all precursors from each multipeptide
            precs =  multipeptides[i].getAllPeptides()
            for prec in precs:
                # Add precursor ID if it is from best_run.
                if prec.getRunId() == run_id:
                    precursor_ids.add(prec.get_id())
        # Get a sorted list
        precursor_ids = sorted(precursor_ids)
        # Assign best_run as reference run for each precursor_id
        reference_run = dict.fromkeys(precursor_ids , self.best_run)
        return reference_run

    # Get a reference run for each precursor
    def _get_precursor_reference_run(self, multipeptides):
        """
        Returns a dectionary with precursor id as key and run as value.
        Precursors may have different reference runs based on FDR score of associated best peak group.
        """
        reference_run = {}
        for i in range(len(multipeptides)):
            # Get precursor groups from each multipeptide
            prec_groups =  multipeptides[i].getPrecursorGroups()
            for prec_group in prec_groups:
                # Get precursors from a precursor group
                precs = prec_group.getAllPrecursors()
                max_fdr = 1.0
                for prec in precs:
                    # Get all precursor ids
                    prec_id = prec.get_id()
                    # Get best FDR value for the precursor
                    cur_fdr = prec.get_best_peakgroup().get_fdr_score()
                    if cur_fdr <= self.alignment_fdr_threshold and cur_fdr < max_fdr:
                        max_fdr = cur_fdr
                        # Make Run of the current precursor as the reference run
                        reference_run[prec_id] = prec.getRun()
                if prec_id not in reference_run:
                    # No peak group has FDR lower than alignment_fdr_threshold
                    reference_run[prec_id] = None
        return reference_run

    # Get a reference run for each precursor-group
    def _get_multipeptide_reference_run(self, multipeptides):
        """
        Returns a dectionary with precursor id as key and run as value.
        Precursors may have different reference runs based on FDR score of associated multipeptide's best peak group.
        """
        reference_run = {}
        for i in range(len(multipeptides)):
            # Get precursor groups from each multipeptide
            prec_groups =  multipeptides[i].getPrecursorGroups()
            max_fdr = 1.0
            refRun = None
            for prec_group in prec_groups:
                precs = prec_group.getAllPrecursors()
                cur_fdr = prec_group.getOverallBestPeakgroup().get_fdr_score()
                if cur_fdr <= self.alignment_fdr_threshold and cur_fdr < max_fdr:
                    max_fdr = cur_fdr
                    refRun = prec_group.run_
                for prec in precs:
                    # Get all precursor ids
                    prec_id = prec.get_id()
                    # Make Run of the current precursor as the reference run
                    reference_run[prec_id] = refRun
        return reference_run
