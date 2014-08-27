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

from msproteomicstoolslib.util.assertions import pre_condition, post_condition, class_invariant
from msproteomicstoolslib.data_structures.PeakGroup import MinimalPeakGroup

class PrecursorGroup():
    """A set of precursors that are isotopically modified versions of each other.

    A collection of precursors that are isotopically modified versions of the
    same underlying peptide sequence. Generally these are heavy/light forms.
    """

    __slots__ = ["peptide_group_label_", "run_", "precursors_"]

    def __init__(self, peptide_group_label, run):
        self.peptide_group_label_ = peptide_group_label  
        self.run_ = run
        self.precursors_ = []

    def __str__(self):
        return "PrecursorGroup %s" % (self.getPeptideGroupLabel())

    def __iter__(self):
        for precursor in self.precursors_:
            yield precursor

    def __classInvariant__(self):
        # for precursor in self.precursors_: print precursor.sequence
        if len(self.precursors_) > 0:
            # All precursor sequences should all be equal to the first sequence
            assert(all( [precursor.sequence == self.precursors_[0].sequence for precursor in self.precursors_] )) 
        return True

    @class_invariant(__classInvariant__)
    def getPeptideGroupLabel(self):
        return self.peptide_group_label_
  
    @class_invariant(__classInvariant__)
    def addPrecursor(self, precursor):
        precursor.precursor_group = self
        self.precursors_.append(precursor)

    @class_invariant(__classInvariant__)
    def getPrecursor(self, curr_id):
        for precursor in self:
            if precursor.get_id() == curr_id:
                return precursor
        return None

    @class_invariant(__classInvariant__)
    def getAllPrecursors(self):
        return list(self)

    @class_invariant(__classInvariant__)
    def getAllPeakgroups(self):
        for pr in self.precursors_:
            for pg in pr.get_all_peakgroups():
                yield pg

    @class_invariant(__classInvariant__)
    def getOverallBestPeakgroup(self):
        allpg = list(self.getAllPeakgroups())
        if len(allpg) == 0:
            return None

        minscore = min([pg.get_fdr_score() for pg in allpg])
        return [pg for pg in allpg if pg.get_fdr_score() <= minscore][0]

