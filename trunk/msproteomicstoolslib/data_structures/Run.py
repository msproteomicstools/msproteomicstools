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

from msproteomicstoolslib.data_structures.PrecursorGroup import PrecursorGroup

class Run():
    """
    A run contains references to identified precursor groups and precursors. 
    
    The run stores a reference to precursor groups (heavy/light pairs) identified in the run.
    It has a unique id and stores the headers from the csv

    A run has the following attributes: 
        - an identifier that is unique to this run
        - a filename where it originally came from
        - a dictionary of precursor groups which are accessible through the following functions
          - getPrecursorGroup 
          - hasPrecursor
          - getPrecursor
          - addPrecursor
    """

    def __init__(self, header, header_dict, runid, orig_input_filename=None, filename=None, aligned_filename=None):
        self.header = header
        self.header_dict = header_dict
        self.runid = runid
        self.orig_filename = orig_input_filename # the original input filename
        self.openswath_filename = filename # the original OpenSWATH filename
        self.aligned_filename = aligned_filename # the aligned filename
        self.all_precursor_groups_ = {}
  
    def __str__(self):
        return "Run %s" % (self.get_id())

    def get_id(self):
        return self.runid

    def get_openswath_filename(self):
        return self.openswath_filename

    def get_aligned_filename(self):
        return self.aligned_filename
  
    def get_best_peaks(self):
        result = []
        for k, precursor_group in self.all_precursor_groups_.iteritems():
            for peptide in precursor_group:
                result.append(peptide.get_best_peakgroup())
        return result
  
    def get_best_peaks_with_cutoff(self, cutoff):
        return [p for p in self.get_best_peaks() if p.get_fdr_score() < cutoff]
  
    def getPrecursorGroup(self, curr_id):
        try:
          return self.all_precursor_groups_[curr_id]
        except KeyError:
          # this run has no peakgroup for that peptide
          return None

    def hasPrecursor(self, peptide_group_label, trgr_id):
        return self.all_precursor_groups_.has_key(peptide_group_label) and \
                not self.getPrecursorGroup(peptide_group_label).getPrecursor(trgr_id) is None

    def getPrecursor(self, peptide_group_label, trgr_id):
        if self.hasPrecursor(peptide_group_label, trgr_id):
            return self.getPrecursorGroup(peptide_group_label).getPrecursor(trgr_id)
        return None

    def addPrecursor(self, precursor, peptide_group_label):
        if self.all_precursor_groups_.has_key(peptide_group_label):
            self.getPrecursorGroup(peptide_group_label).addPrecursor(precursor)
        else:
            prec_gr = PrecursorGroup(peptide_group_label, self)
            prec_gr.addPrecursor(precursor)
            self.all_precursor_groups_[peptide_group_label] = prec_gr

    def __iter__(self):
        for precursor in self.all_precursor_groups_.values():
            yield precursor

