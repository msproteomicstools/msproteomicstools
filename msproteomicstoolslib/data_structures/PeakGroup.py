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

class PeakGroupBase(object):

    __slots__ = ["fdr_score", "normalized_retentiontime", "id_", "intensity_", "cluster_id_"]

    def __init__(self):
        self.fdr_score = None
        self.normalized_retentiontime = None
        self.id_ = None
        self.intensity_ = None
        self.cluster_id_ = -1
  
    def __str__(self):
        return "PeakGroup %s at %s s in %s with score %s (cluster %s)" % (self.get_feature_id(), self.get_normalized_retentiontime(), self.peptide.run, self.get_fdr_score(), self.get_cluster_id())

    def get_value(self, value):
        raise Exception("Needs implementation")

    def set_value(self, key, value):
        raise Exception("Needs implementation")

    def set_fdr_score(self, fdr_score):
        self.fdr_score = fdr_score

    def get_fdr_score(self):
        return self.fdr_score

    def set_normalized_retentiontime(self, normalized_retentiontime):
        self.normalized_retentiontime = float(normalized_retentiontime)

    def get_normalized_retentiontime(self):
        return self.normalized_retentiontime

    def set_feature_id(self, id_):
        self.id_ = id_

    def get_feature_id(self):
        return self.id_

    def set_intensity(self, intensity):
      self.intensity_ = intensity

    def get_intensity(self):
        return self.intensity_

    # Selected
    def is_selected(self):
        return self.cluster_id_ == 1

    def select_this_peakgroup(self):
        self.cluster_id_ = 1

    def get_cluster_id(self):
        return self.cluster_id_

class MinimalPeakGroup(PeakGroupBase):
    """
    A single peakgroup that is defined by a retention time in a chromatogram
    of multiple transitions. Additionally it has an fdr_score and it has an
    aligned RT (e.g. retention time in normalized space).
    A peakgroup can be selected for quantification or not (this is stored as
    having cluster_id == 1).
    
    Note that for performance reasons, the peakgroups are created on-the-fly
    and not stored as objects but rather as tuples in "Peptide".

    Each peak group has a unique id, a score (fdr score usually), a retention
    time as well as a back-reference to the precursor that generated the
    peakgroup.
    In this case, the peak group can also be assigned a cluster id (where the
    cluster 1 is special as the one we will use for quantification).
    """

    def __init__(self, unique_id, fdr_score, assay_rt, selected, cluster_id, peptide, intensity=None, dscore=None):
      super(MinimalPeakGroup, self).__init__()
      self.id_ = unique_id
      self.fdr_score = fdr_score
      self.normalized_retentiontime = assay_rt 
      self.cluster_id_ = cluster_id
      self.peptide = peptide
      self.intensity_ = intensity
      self.dscore_ = dscore
  
    ## Print
    def print_out(self):
        # return self.run.get_id() + "/" + self.get_id() + " " + str(self.get_fdr_score()) + " " + str(self.get_normalized_retentiontime()) + " " + str(self.get_value("RT")) + " " + str(self.get_value("rt_score")) # rt_score = delta iRT
        return self.peptide.run.get_id() + "/" + self.get_feature_id() + " score:" + str(self.get_fdr_score()) + " RT:" + str(self.get_normalized_retentiontime()) # + " " + str(self.get_value("RT")) + " " + str(self.get_value("rt_score")) # rt_score = delta iRT

    # Do not allow setting of any parameters (since data is not stored here)
    def set_fdr_score(self, fdr_score):
        raise Exception("Cannot set in immutable object")

    def set_normalized_retentiontime(self, normalized_retentiontime):
        raise Exception("Cannot set in immutable object")

    def set_feature_id(self, id_):
        raise Exception("Cannot set in immutable object")

    def set_intensity(self, intensity):
        raise Exception("Cannot set in immutable object")

    def get_dscore(self):
        return self.dscore_

    ## Select / De-select peakgroup
    def select_this_peakgroup(self):
        self.peptide.select_pg(self.get_feature_id())

    ## Select / De-select peakgroup
    def setClusterID(self, id_):
        self.cluster_id_ = id_
        self.peptide.setClusterID(self.get_feature_id(), id_)

    def get_cluster_id(self):
        return self.cluster_id_

class GuiPeakGroup(PeakGroupBase):
    """
    A single peakgroup that is defined by a retention time in a chromatogram
    of multiple transitions.
    """
    def __init__(self, fdr_score, intensity, leftWidth, rightWidth, peptide):
      super(PeakGroupBase, self).__init__()
      self.fdr_score = fdr_score
      self.intensity_ = intensity
      self.leftWidth_ = leftWidth
      self.rightWidth_ = rightWidth
      self.peptide = peptide
  
    def get_value(self, value):
        if value == "m_score":
            return self.fdr_score
        elif value == "Intensity":
            return self.intensity_
        elif value == "rightWidth":
            return self.rightWidth_
        elif value == "leftWidth":
            return self.leftWidth_
        elif value == "FullPeptideName":
            return self.peptide.sequence
        elif value == "Charge":
            return self.charge
        else:
            raise Exception("Do not have value " + value)

class GeneralPeakGroup(PeakGroupBase):

    __slots__ = ["row", "run", "peptide"]

    def __init__(self, row, run, peptide):
      super(GeneralPeakGroup, self).__init__()
      self.row = row
      self.run = run
      self.peptide = peptide

    def get_value(self, value):
        return self.row[self.run.header_dict[value]]

    def set_value(self, key, value):
        if value is None:
            value = "NA"
        self.row[self.run.header_dict[key]] = value

    def get_dscore(self):
        return self.get_value("d_score")

    def print_out(self):
        return self.peptide.run.get_id() + "/" + self.get_feature_id() + " score:" + str(self.get_fdr_score()) + " norm_RT:" + str(self.get_normalized_retentiontime()) + " RT:" + str(self.get_value("RT")) + " Int : " + str(self.get_value("Intensity"))
  
    def setClusterID(self, clid):
        self.cluster_id_ = clid

