#!/usr/bin/python
# -*- coding: utf-8 -*-
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

import os

from SqlDataAccess import SqlDataAccess
from FormatHelper import FormatHelper

class SqlSwathRun():
    """Data Model for a single sqMass file.

    TODO: each file may contain multiple runs!

    Attributes:
        runid: Current run id

    Private Attributes:
       - _run:        A :class:`.SqlDataAccess` object
       - _filename:   Original filename
       - _basename:   Original filename basename
       - _precursor_mapping:   Dictionary { FullPrecursorName : [transition_id, transition_id] }
       - _sequences_mapping:   Dictionary { StrippedSequence : [FullPrecursorName, FullPrecursorName]}

    """

    def __init__(self, runid, filename, load_in_memory=False, precursor_mapping = None, sequences_mapping = None, protein_mapping = {}):
        print "runid ", runid
        if runid is not None:
            assert len(runid) == 1

        self.runid = runid[0] # TODO 
        self._filename = filename
        self._basename = os.path.basename(filename)

        self._range_mapping = {}
        self._score_mapping = {}
        self._intensity_mapping = {}
        self._assay_mapping = {}

        self._run = SqlDataAccess(filename)

        # Map which holds the relationship between a precursor and its
        # corresponding chromatogram id (transition id).
        self._precursor_mapping = {}
        # Map which holds the relationship between a sequence and its precursors
        self._sequences_mapping = {}
        # Map which holds the relationship between a protein and its sequences 
        self._protein_mapping = protein_mapping
        # Map which holds the relationship between a chromatogram and its id
        self._id_mapping = {}

        self._in_memory = False

        if not precursor_mapping is None and not len(precursor_mapping) == 0 \
          and not sequences_mapping is None and not len(sequences_mapping) == 0:
            self._precursor_mapping = precursor_mapping
            self._sequences_mapping = sequences_mapping
            self._id_mapping = dict(self._get_id_mapping())
        else:
            self._group_by_precursor()
            self._group_precursors_by_sequence()

    def get_transitions_for_precursor_display(self, precursor):

        if not self._precursor_mapping.has_key(str(precursor)):
            return [ "NA" ]

        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            transitions.append(chrom_id)
        return transitions

    def get_range_data(self, precursor):
        return self._range_mapping.get(precursor, [ [0,0] ])

    def get_assay_data(self, precursor):
        r = self._assay_mapping.get(precursor, None)

        if r is not None and len(r) > 0:
            return r[0]
        return None

    def get_score_data(self, precursor):
        r = self._score_mapping.get(precursor, None)

        if r is not None and len(r) > 0:
            return r[0]
        return None

    def get_intensity_data(self, precursor):
        r = self._intensity_mapping.get(precursor, None)

        if r is not None and len(r) > 0:
            return r[0]
        return None

    #
    ## Initialization
    #

    def _get_id_mapping(self):
        """
        Map SQL ids to transition group identifiers
        """

        import sqlite3
        conn = sqlite3.connect(self._filename)
        c = conn.cursor()

        id_mapping = [row for row in c.execute("SELECT NATIVE_ID, ID FROM CHROMATOGRAM" )]

        return id_mapping

    def _group_by_precursor(self):
        """
        Populate the mapping between precursors and the chromatogram ids.

        The precursor is of type 'PEPT[xx]IDE/3' with an optional (DECOY) tag
        in front. Different modifications will generate different precursors,
        but not different charge states.
        """

        id_mapping = self._get_id_mapping()
        # TODO: could also be done in SQL
        self._id_mapping = dict(id_mapping)

        openswath_format = self._has_openswath_format([m[0] for m in id_mapping])

        result = {}
        if not openswath_format:
            raise Exception("Could not parse chromatogram ids ... ")

        f = FormatHelper()
        for key, myid in id_mapping:
            trgr_nr = f._compute_transitiongroup_from_key(key)
            tmp = self._precursor_mapping.get(trgr_nr, [])
            tmp.append(key)
            self._precursor_mapping[trgr_nr] = tmp

    def _has_openswath_format(self, keys):
        f = FormatHelper()
        return all([ f._has_openswath_format(k) for k in keys] )

    def _group_precursors_by_sequence(self):
        """Group together precursors with the same charge state"""
        self._sequences_mapping = {}
        for precursor in self._precursor_mapping.keys():
            seq = precursor.split("/")[0]
            tmp = self._sequences_mapping.get(seq, [])
            tmp.append(precursor)
            self._sequences_mapping[seq] = tmp

    #
    ## Getters (data) -> see ChromatogramTransition.getData
    #

    def getTransitionCount(self):
        """
        Get total number of transitions 
        """
        return sum([ len(v) for v in self._precursor_mapping.values()])

    def get_data_for_precursor(self, precursor):
        """Retrieve raw data for a specific precursor - data will be as list of
        pairs (timearray, intensityarray)"""

        if not self._precursor_mapping.has_key(str(precursor)):
            return [ [ [0], [0] ] ]

        transitions = []
        sql_ids = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            sql_id = self._id_mapping[ chrom_id ]
            sql_ids.append(sql_id)
        transitions = self._run.getDataForChromatograms(sql_ids)
        return transitions

    def get_data_for_transition(self, transition_id):
        """
        Retrieve raw data for a specific transition
        """

        if transition_id in self._id_mapping:
            return [self._run.getDataForChromatogram(self._id_mapping[transition_id])]
        else:
            print "Warning: Found chromatogram identifier '%s' that does not map to any chromatogram in the data." % transition_id
            print "Please check your input data"

    def get_id(self):
        return self._basename

    #
    ## Getters (info)
    #
    def get_transitions_for_precursor(self, precursor):
        """
        Return the transition names for a specific precursor
        """
        return self._precursor_mapping.get(str(precursor), [])

    def get_sequence_for_protein(self, protein):
        return self._protein_mapping.get(protein, [])

    def get_precursors_for_sequence(self, sequence):
        """
        Get all precursors mapping to one stripped sequence
        """
        return self._sequences_mapping.get(sequence, [])

    def get_all_precursor_ids(self):
        """
        Get all precursor ids (full sequence + charge)
        """
        return self._precursor_mapping.keys()

    def get_all_peptide_sequences(self):
        """
        Get all (stripped) sequences
        """
        return self._sequences_mapping.keys()

    def get_all_proteins(self):
        return self._protein_mapping

    # 
    ## Data manipulation
    #
    def remove_precursors(self, toremove):
        """ Remove a set of precursors from the run (this can be done to filter
        down the list of precursors to display).
        """
        for key in toremove:
            self._precursor_mapping.pop(key, None)
        self._group_precursors_by_sequence()

        # Re-initialize self to produce correct mapping
        ## self._initialize()

    def add_peakgroup_data(self, precursor_id, leftWidth, rightWidth, fdrscore, intensity, assay_rt):

        tmp = self._range_mapping.get(precursor_id, [])
        tmp.append( [leftWidth, rightWidth ] )
        self._range_mapping[precursor_id] = tmp

        tmp = self._score_mapping.get(precursor_id, [])
        tmp.append(fdrscore)
        self._score_mapping[precursor_id] = tmp

        tmp = self._intensity_mapping.get(precursor_id, [])
        tmp.append(intensity)
        self._intensity_mapping[precursor_id] = tmp

        tmp = self._assay_mapping.get(precursor_id, [])
        tmp.append(assay_rt)
        self._assay_mapping[precursor_id] = tmp

