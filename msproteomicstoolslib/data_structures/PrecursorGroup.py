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

class PrecursorGroup(object):
    """A set of precursors that are isotopically modified versions or different
    charge states of each other.

    A collection of precursors that are isotopically modified versions or
    different charge states of the same underlying peptide sequence. Generally
    these are heavy/light forms. This class groups these Precursors together.

    Attributes:
        - self.peptide_group_label_: Identifier or precursor group 
        - self.run_: Reference to the :class:`.Run` where this PrecursorGroup is from
        - self.precursors_: List of actual precursors
    """

    __slots__ = ["peptide_group_label_", "run_", "precursors_"]

    def __init__(self, peptide_group_label, run):
        self.peptide_group_label_ = peptide_group_label  
        self.run_ = run
        self.precursors_ = []

    def __str__(self):
        return "PrecursorGroup %s" % (self.getPeptideGroupLabel())

    def __lt__(self, other):

        if self.run_.get_id() == other.run_.get_id():
            return self.getPeptideGroupLabel() > other.getPeptideGroupLabel()
        else:
            return self.run_.get_id() > other.run_.get_id()

    def __iter__(self):
        for precursor in self.precursors_:
            yield precursor

    def __classInvariant__(self):
        # for precursor in self.precursors_: print precursor.sequence
        if len(self.precursors_) > 0:
            # All precursor sequences should all be equal to the first sequence
            assert(all( [precursor.getSequence() == self.precursors_[0].getSequence() for precursor in self.precursors_] )) 
        return True

    @class_invariant(__classInvariant__)
    def getPeptideGroupLabel(self):
        """
        getPeptideGroupLabel(self)
        Get peptide group label
        """
        return self.peptide_group_label_
  
    @class_invariant(__classInvariant__)
    def addPrecursor(self, precursor):
        """
        addPrecursor(self, precursor)
        Add precursor to peptide group
        """
        precursor.set_precursor_group( self )
        self.precursors_.append(precursor)

    @class_invariant(__classInvariant__)
    def getPrecursor(self, curr_id):
        """
        getPrecursor(self, curr_id)
        Get the precursor for the given transition group id
        """
        for precursor in self:
            if precursor.get_id() == curr_id:
                return precursor
        return None

    @class_invariant(__classInvariant__)
    def getAllPrecursors(self):
        """
        getAllPrecursors(self)
        Return a list of all precursors in this precursor group
        """
        return list(self)

    @class_invariant(__classInvariant__)
    def getAllPeakgroups(self):
        """
        getAllPeakgroups(self)
        Generator of all peakgroups attached to the precursors in this group
        """
        for pr in self.precursors_:
            for pg in pr.get_all_peakgroups():
                yield pg

    @class_invariant(__classInvariant__)
    def getOverallBestPeakgroup(self):
        """
        getOverallBestPeakgroup(self)
        Get the best peakgroup (by fdr score) of all precursors contained in this precursor group
        """
        allpg = list(self.getAllPeakgroups())
        if len(allpg) == 0:
            return None

        minscore = min([pg.get_fdr_score() for pg in allpg])
        return [pg for pg in allpg if pg.get_fdr_score() <= minscore][0]

    def get_decoy(self):
        """
        Whether the current peptide is a decoy or not

        Returns:
            decoy(bool): Whether the peptide is decoy or not
        """
        if len(self.precursors_) == 0:
            return False

        return self.precursors_[0].get_decoy()

