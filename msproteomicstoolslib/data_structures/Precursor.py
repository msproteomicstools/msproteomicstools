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

from __future__ import print_function
from msproteomicstoolslib.data_structures.PeakGroup import MinimalPeakGroup

class PrecursorBase(object):

    __slots__ = ["_decoy"]

    def __init__(self, this_id, run):
        raise NotImplemented

    def get_id(self):
        return self.id 
  
    def get_decoy(self):
        return self._decoy

    def set_decoy(self, decoy):
        if decoy in ["FALSE", "False", "0"]:
            self._decoy = False
        elif decoy in ["TRUE", "True", "1"]:
            self._decoy = True
        else:
            raise Exception("Unknown decoy classifier '%s', please check your input data!" % decoy)
  
    # store information about the peakgroup - tuples (e.g. whether they are selected)
    def select_pg(self, this_id):
        raise NotImplemented

    def unselect_pg(self, id):
        raise NotImplemented

    def get_best_peakgroup(self):
        raise NotImplemented

    def get_selected_peakgroup(self):
        raise NotImplemented

    def get_all_peakgroups(self):
        raise NotImplemented
  
    def find_closest_in_iRT(self, delta_assay_rt):
        raise NotImplemented

class GeneralPrecursor(PrecursorBase):
    """ A set of peakgroups that belong to the same precursor in a single run.

    == Implementation details ==
    
    This is a plain implementation where all peakgroup objects are stored in a
    simple list, this is not very efficient since many objects need to be
    created which in Python takes a lot of memory.
    """

    __slots__ = ["id", "run", "peakgroups", "protein_name", "sequence", "precursor_group"]

    def __init__(self, this_id, run):
        self.id = this_id
        self.peakgroups = []
        self.run = run
        self._decoy = False
        self.protein_name = ""
        self.sequence = ""
        self.precursor_group = None
  
    def __str__(self):
        return "%s (run %s)" % (self.id, self.run)

    def add_peakgroup(self, peakgroup):
        self.peakgroups.append(peakgroup)
  
    def getRun(self):
      return self.run
  
    def set_precursor_group(self, p):
        self.precursor_group = p

    def getSequence(self):
        return self.sequence

    def setSequence(self, s):
        self.sequence = s
  
    def setProteinName(self, p):
        self.protein_name = p 

    def getProteinName(self):
        return self.protein_name

    def get_run_id(self):
      return self.run.get_id()

    def getRunId(self):
      return self.run.get_id()
  
    def append(self, transitiongroup):
        assert self.id == transitiongroup.get_id()
        self.peakgroups.append(transitiongroup)
  
    def get_best_peakgroup(self):
        """ Return the best peakgroup according to fdr score
        """
        if len(self.peakgroups) == 0: return None
        best_score = self.peakgroups[0].get_fdr_score()
        result = self.peakgroups[0]
        for peakgroup in self.peakgroups:
            # TODO this is not well done !!! 
            if best_score is None or \
               peakgroup.get_fdr_score() is None or \
            peakgroup.get_fdr_score() <= best_score:
                # print "better: ", peakgroup.get_fdr_score(),  best_score
                best_score = peakgroup.get_fdr_score()
                result = peakgroup
        return result
  
    def get_selected_peakgroup(self):
        # return the selected peakgroup of this precursor, we can only select 1 or
        # zero groups per chromatogram!
        selected = [peakgroup for peakgroup in self.peakgroups if peakgroup.is_selected()]
        assert len(selected) < 2
        if len(selected) == 1:
          return selected[0]
        else: 
            return None

    def get_all_peakgroups(self):
        return self.peakgroups
  
    def find_closest_in_iRT(self, delta_assay_rt):
      return min(self.peakgroups, key=lambda x: abs(float(x.get_normalized_retentiontime()) - float(delta_assay_rt)))

class Precursor(PrecursorBase):
    """ A set of peakgroups that belong to the same precursor in a single run.

    Each precursor has a backreference to its precursor group (heavy/light
    pair) it belongs to, the run it belongs to as well as its amino acid sequence.
    Furthermore, a unique id for the precursor and the protein name are stored.

    A precursor can return its best transition group, the selected peakgroup,
    or can return the transition group that is closest to a given iRT time.
    Its id is the transition_group_id (e.g. the id of the chromatogram)

    The "selected" peakgroup is represented by the peakgroup that belongs to
    cluster number 1 (cluster_id == 1) which in this case is "special".

    == Implementation details ==
    
    For memory reasons, we store all information about the peakgroup in a
    tuple (invariable). This tuple contains a unique feature id, a score and
    a retention time. Additionally, we also store, in which cluster the
    peakgroup belongs (if the user sets this).

    A peakgroup has the following attributes: 
        - an identifier that is unique among all other precursors 
        - a set of peakgroups 
        - a back-reference to the run it belongs to
    """

    __slots__ = ["id", "run", "peakgroups_", "cluster_ids_", "protein_name", "sequence", "precursor_group"]

    def __init__(self, this_id, run):
        self.id = this_id  
        self.run = run
        self.peakgroups_ = []
        self.cluster_ids_ = []
        self._decoy = False
        self.protein_name = ""
        self.sequence = ""
        self.precursor_group = None
  
    def getRun(self):
      return self.run

    def __str__(self):
        return "%s (run %s)" % (self.id, self.run)

    def add_peakgroup_tpl(self, pg_tuple, tpl_id, cluster_id=-1):
        """Adds a peakgroup to this precursor.

        The peakgroup should be a tuple of length 4 with the following components:
            0. id
            1. quality score (FDR)
            2. retention time (normalized)
            3. intensity
            (4. d_score optional)
        """
        # Check that the peak group is added to the correct precursor
        if self.id != tpl_id:
            raise Exception("Cannot add a tuple to this precursor with a different id")

        if len(pg_tuple) == 4:
            pg_tuple = pg_tuple + (None,)

        assert len(pg_tuple) == 5
        self.peakgroups_.append(pg_tuple)
        self.cluster_ids_.append(cluster_id)

    def get_id(self):
        return self.id 

    def set_precursor_group(self, p):
        self.precursor_group = p

    def getSequence(self):
        return self.sequence

    def setSequence(self, s):
        self.sequence = s
  
    def setProteinName(self, p):
        self.protein_name = p 

    def getProteinName(self):
        return self.protein_name

    def getPrecursorGroup(self):
        return self.precursor_group 

    def getRunId(self):
      return self.run.get_id()

    def get_run_id(self):
      return self.run.get_id()
    
    # 
    # Peakgroup cluster membership
    # 
    def select_pg(self, this_id):
        self.setClusterID(this_id, 1)

    def unselect_pg(self, this_id):
        self.setClusterID(this_id, -1)

    def setClusterID(self, this_id, cl_id):
        pg_id = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == this_id]
        assert len(pg_id) == 1
        assert cl_id == -1 or len([0 for i in self.cluster_ids_ if i == cl_id] ) == 0
        self.cluster_ids_[pg_id[0]] = cl_id

    def unselect_all(self):
        for i in range(len(self.cluster_ids_)) : 
            self.cluster_ids_[i] = -1

    # 
    # Peakgroup selection
    # 
    def get_best_peakgroup(self):
        if len(self.peakgroups_) == 0:
            return None

        best_score = self.peakgroups_[0][1]
        result = self.peakgroups_[0]
        for peakgroup in self.peakgroups_:
            if peakgroup[1] <= best_score:
                best_score = peakgroup[1]
                result = peakgroup
        index = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == result[0]][0]
        return MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])

    def _fixSelectedPGError(self, fixMethod="Exception"):
      selected = [i for i,pg in enumerate(self.cluster_ids_) if pg == 1]
      if len(selected) > 1:
          print("Potential error detected in %s:\nWe have more than one selected peakgroup found. Starting error handling by using method '%s'." % (self, fixMethod))
          best_score = self.peakgroups_[0][1]
          best_pg = 0
          for s in selected:
              if best_score > self.peakgroups_[s][1]:
                  best_score = self.peakgroups_[s][1]
                  best_pg = s

          if fixMethod == "Exception":
              raise Exception("More than one selected peakgroup found in %s " % self )
          elif fixMethod == "BestScore":
              # Deselect all, then select the one with the best score...
              for s in selected:
                  self.cluster_ids_[s] = -1
              self.cluster_ids_[best_pg] = 1


    def get_selected_peakgroup(self):
      # return the selected peakgroup of this precursor, we can only select 1 or
      # zero groups per chromatogram!
      selected = [i for i,pg in enumerate(self.cluster_ids_) if pg == 1]
      assert len(selected) < 2
      if len(selected) == 1:
        index = selected[0]
        result = self.peakgroups_[index]
        return MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])
      else: 
          return None

    def getClusteredPeakgroups(self):
      selected = [i for i,pg in enumerate(self.cluster_ids_) if pg != -1]
      for index in selected:
        result = self.peakgroups_[index]
        yield MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])

    def get_all_peakgroups(self):
        for index, result in enumerate(self.peakgroups_):
            yield MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])

    def getAllPeakgroups(self):
        return self.get_all_peakgroups()
  
    def find_closest_in_iRT(self, delta_assay_rt):
      result = min(self.peakgroups_, key=lambda x: abs(float(x[2]) - float(delta_assay_rt)))
      index = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == result[0]][0]
      return MinimalPeakGroup(result[0], result[1], result[2], self.cluster_ids_[index] == 1, self.cluster_ids_[index], self, result[3], result[4])

