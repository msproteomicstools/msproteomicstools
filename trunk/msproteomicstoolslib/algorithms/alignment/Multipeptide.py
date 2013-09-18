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

class Multipeptide():
    """ A collection of the same precursors (chromatograms) across multiple runs.

    It contains individual precursors that can be accessed by their run id.
    """
  
    def __init__(self):
        self._peptides = {}
        self._has_null = False
        self._nr_runs = -1

    def __str__(self):
        return "Precursors of %s runs, identified by %s." % (len(self._nr_runs), self.get_peptides()[0].id)
  
    # 
    ## Getters  / Setters
    # 
    def set_nr_runs(self, v):
        self._nr_runs = v

    def get_nr_runs(self):
        return self._nr_runs

    def has_peptide(self, runid):
        return self._peptides.has_key(runid)

    def get_peptide(self, runid):
        return self._peptides[runid]

    def get_peptides(self):
      return self._peptides.values()

    def get_id(self):
      if len(self.get_peptides()) == 0: return None
      return self.get_peptides()[0].get_id()

    def more_than_fraction_selected(self, fraction):
      assert self._nr_runs >= 0
      # returns true if more than fraction of the peakgroups are selected
      if len( self.get_selected_peakgroups() )*1.0 / self._nr_runs < fraction:
          return False
      return True

    def get_decoy(self):
        if len(self.get_peptides()) == 0: return False
        return self.get_peptides()[0].get_decoy() 

    def has_null_peptides(self):
      return self._has_null

    def insert(self, runid, peptide):
      assert not self.has_peptide(runid)

      if peptide is None: 
          self._has_null = True 
          return
      self._peptides[runid] = peptide
    
    def get_selected_peakgroups(self):
      return [p.get_selected_peakgroup() for p in self.get_peptides() if p.get_selected_peakgroup() is not None]

    def find_best_peptide_pg(self):
      # Find best peakgroup across all peptides
      best_fdr = 1.0
      for p in self.get_peptides():
        if(p.get_best_peakgroup().get_fdr_score() < best_fdr): 
            result = p.get_best_peakgroup()
            best_fdr = p.get_best_peakgroup().get_fdr_score() 
      return result
  
    # 
    ## Methods
    #

    def detect_outliers(self):
        from msproteomicstoolslib.math.chauvenet import chauvenet
        import numpy
        # Uses chauvenet's criterion for outlier detection to find peptides
        # whose retention time is different from the rest.
        rts = [float(p.get_selected_peakgroup().get_normalized_retentiontime()) for p in self.get_peptides() if p.get_selected_peakgroup() is not None]
        runids = numpy.array([p.get_run_id() for p in self.get_peptides() if p.get_selected_peakgroup() is not None])
        if len(rts) == 1: return []
        outliers = chauvenet(numpy.array(rts),numpy.array(rts))
        return runids[~outliers]

    # 
    ## Boolean questions
    #

    def all_above_cutoff(self, cutoff):
      assert self._nr_runs >= 0
      if len(self.get_peptides())< self._nr_runs:
          return False

      for p in self.get_peptides():
        if p.get_best_peakgroup().get_fdr_score() > cutoff: 
            return False
      return True
  
    def all_selected(self):
      assert self._nr_runs >= 0
      if len(self.get_peptides())< self._nr_runs:
          return False

      for p in self.get_peptides():
          if p.get_selected_peakgroup() is None: return False
      return True

