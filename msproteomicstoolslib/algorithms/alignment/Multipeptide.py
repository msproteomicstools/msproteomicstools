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

from msproteomicstoolslib.data_structures.PrecursorGroup import PrecursorGroup

class Multipeptide(object):
    """ A collection of the same precursors (chromatograms) across multiple runs.

    It contains individual precursors that can be accessed by their run id.
    """
  
    def __init__(self):
        self._peptides = {}
        self._has_null = False
        self._nr_runs = -1

    def __str__(self):
        if len(self.getPrecursorGroups()) > 0:
            return "Precursors of %s runs, identified by %s." % (
                self._nr_runs, list(self.getPrecursorGroups())[0].getPeptideGroupLabel())
        else:
            return "Empty set of precursors."
  
    # 
    ## Getters  / Setters
    # 
    def set_nr_runs(self, v):
        self._nr_runs = v

    def get_nr_runs(self):
        return self._nr_runs

    def has_peptide(self, runid):
        raise Exception("This doesnt do what you want")
        return runid in self._peptides

    def get_peptide(self, runid):
        raise Exception("This doesnt do what you want")
        return self._peptides[runid]

    def get_peptides(self):
        raise Exception("This doesnt do what you want")
        return self._peptides.values()

    def getPrecursorGroup(self, runid):
        """
        Get precursor group for the given run

        Parameters
        ----------
        :param str runid: Run id of the group

        :rtype: :class:`.PrecursorGroup`: Precursor group from the corresponding run
        """
        return self._peptides[runid]

    def getPrecursorGroups(self):
        """
        Get all precursor groups

        :rtype: list(:class:`.PrecursorGroup`): All Precursor group from the corresponding run
        """
        return sorted(self._peptides.values())

    def hasPrecursorGroup(self, runid):
        """
        Checks whether a given run has a precursor group

        :param str runid: Run id to check
        :rtype: bool: Whether the given run has a precursor group
        """
        return runid in self._peptides

    def getAllPeptides(self):
      return [p for prgr in self.getPrecursorGroups() for p in prgr]

    def get_id(self):
      if len(self.getAllPeptides()) == 0:
           return None

      return self.getAllPeptides()[0].get_id()

    def more_than_fraction_selected(self, fraction):
      assert self._nr_runs >= 0
      # returns true if more than fraction of the peakgroups are selected
      if len( self.get_selected_peakgroups() ) *1.0 / self._nr_runs < fraction:
          return False
      return True

    def get_decoy(self):
        """
        Whether the current peptide is a decoy or not

        :rtype: bool: Whether the peptide is decoy or not
        """
        if len(self._peptides) == 0:
            return False

        return next(iter(self._peptides.values())).get_decoy()

    def has_null_peptides(self):
        """
        Whether there are runs in which no peptide was detected (peptide is Null)

        Returns:
            has_null(bool): Whether there are Null peptides in this object (not detected in some runs)
        """
        return self._has_null

    def insert(self, runid, precursor_group):
        """
        Insert a :class:`.PrecursorGroup` into the Multipeptide

        Args:
            runid(str): Run id of the group
            precursor_group(:class:`.PrecursorGroup`): Precursor group to be inserted

        Raises:
            Exception: If self.hasPrecursorGroup(runid) is true
        """

        # Deal with None (store that we have a None precursor)
        if precursor_group is None: 
            self._has_null = True 
            return

        if self.hasPrecursorGroup(runid):
            raise Exception("A precursor for run %s already exists, cannot add another one.")

        self._peptides[runid] = precursor_group

    def get_selected_peakgroups(self):
        """
        Get all peakgroups that were selected across all runs and precursor groups
        """
        return [precursor.get_selected_peakgroup() for prgr in self.getPrecursorGroups() for precursor in prgr if precursor.get_selected_peakgroup() is not None]

    def find_best_peptide_pg(self):
      """
      Find best peakgroup across all peptides
      """
      best_fdr = 1.0
      result = None
      for p in self.getAllPeptides():
        if p.get_best_peakgroup().get_fdr_score() < best_fdr: 
            result = p.get_best_peakgroup()
            best_fdr = p.get_best_peakgroup().get_fdr_score() 
      return result
  
    # 
    ## Boolean questions
    #

    def all_above_cutoff(self, cutoff):
      assert self._nr_runs >= 0
      if len(self.getPrecursorGroups()) < self._nr_runs:
          return False

      for prgr in self.getPrecursorGroups():
          for p in prgr:
              if p.get_best_peakgroup().get_fdr_score() > cutoff: 
                  return False
      return True
  
    def all_selected(self):
        """
        Returns True if all peakgroups are selected
        """
        assert self._nr_runs >= 0
        if len(self.getAllPeptides()) < self._nr_runs:
            return False

        for p in self.getAllPeptides():
            if p.get_selected_peakgroup() is None:
                return False
        return True

