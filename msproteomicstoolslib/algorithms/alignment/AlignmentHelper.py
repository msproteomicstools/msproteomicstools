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

class Multipeptide():
    """
    A collection of the same precursors (chromatograms) across multiple runs.

    It contains individual precursors that can be accessed by their run id.
    """
  
    def __init__(self):
        self._peptides = {}
        self._has_null = False

    def __str__(self):
        return "Precursors of %s runs, identified by %s." % (len(self._peptides), self.get_peptides()[0].id)
  
    # 
    ## Getters  / Setters
    # 

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
      # returns true if more than fraction of the peakgroups are selected
      if len( self.get_selected_peakgroups() )*1.0 / len(self._peptides) < fraction:
          return False
      return True

    def get_decoy(self):
        if len(self.get_peptides()) == 0: return False
        return self.get_peptides()[0].get_decoy() 

    def has_null_peptides(self):
      return self._has_null

    def insert(self, runid, peptide):
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
        # Uses chauvenet's criterion for outlier detection to find peptides
        # whose retention time is different from the rest.
        rts = [float(p.get_selected_peakgroup().get_normalized_retentiontime()) for p in self.get_peptides() if p.get_selected_peakgroup() is not None]
        runids = numpy.array([p.get_selected_peakgroup().get_run_id() for p in self.get_peptides() if p.get_selected_peakgroup() is not None])
        if len(rts) == 1: return []
        outliers = chauvenet(numpy.array(rts),numpy.array(rts))
        return runids[~outliers]

    # 
    ## Boolean questions
    #

    def all_above_cutoff(self, cutoff):
      for p in self.get_peptides():
        if p.get_best_peakgroup().get_fdr_score() > cutoff: 
            return False
      return True
  
    def all_below_cutoff(self, cutoff):
      for p in self.get_peptides():
        if p.get_best_peakgroup().get_fdr_score() < cutoff: return False
      return True

    def all_selected(self):
      for p in self.get_peptides():
          if p.get_selected_peakgroup() is None: return False
      return True

class AlignmentExperiment(object):
    """
    An Alignment Experiment is a container for multiple experimental runs - some of which may contain the same precursors.
    """

    def __init__(self):
        self.runs = []

    def get_all_multipeptides(self, fdr_cutoff, verbose=False):
        # Find all precursors that are above the fdr cutoff in each run and
        # build a union of those precursors. Then search for each of those
        # precursors in all the other runs and build a multipeptide /
        # multiprecursor.
        union_transition_groups = []
        union_proteins = []
        union_target_transition_groups = []
        for i,r in enumerate(self.runs):
            if verbose: 
                stdout.write("\rParsing run %s out of %s" % (i+1, len(self.runs) ))
                stdout.flush()
            union_target_transition_groups.append( [peak.peptide.get_id() for peak in r.get_best_peaks_with_cutoff(fdr_cutoff) if not peak.peptide.get_decoy()] )
            union_transition_groups.append( [peak.peptide.get_id() for peak in r.get_best_peaks_with_cutoff(fdr_cutoff)] )
            union_proteins.append( list(set([peak.peptide.protein_name for peak in r.get_best_peaks_with_cutoff(fdr_cutoff) if not peak.peptide.get_decoy()])) )
        if verbose: stdout.write("\r\r\n") # clean up

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

        if verbose:
            print "==================================="
            print "Finished parsing, number of precursors and peptides per run"
            print "All precursors", [len(s) for s in union_transition_groups], "(union of all runs %s)" % len(self.union_transition_groups_set)
            print "All target precursors", [len(s) for s in union_target_transition_groups], "(union of all runs %s)" % len(union_target_transition_groups_set)
            print "All target proteins", [len(s) for s in union_proteins], "(union of all runs %s)" % len(self.union_proteins_set)
            print "Decoy percentage on precursor level %0.4f%%" % ( (all_prec - target_prec) * 100.0 / all_prec )

        self.estimated_decoy_pcnt =  (all_prec - target_prec) * 100.0 / all_prec 
        if all_prec - target_prec == 0: self.estimated_decoy_pcnt = None

        multipeptides = []
        for peptide_id in self.union_transition_groups_set:
          m = Multipeptide()
          for r in self.runs:
            m.insert(r.get_id(), r.get_peptide(peptide_id))
          multipeptides.append(m)
        return multipeptides


