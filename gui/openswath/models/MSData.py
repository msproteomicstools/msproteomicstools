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


class RunDataModel():

    def __init__(self, run, filename):
        import os
        self._run = run
        self._filename = filename
        self._basename = os.path.basename(filename)
        self._precursor_mapping = {}
        self._sequences_mapping = {}
        self._range_mapping = {}
        self._score_mapping = {}
        self._intensity_mapping = {}

        self._group_by_precursor()
        self._group_precursors_by_sequence()

    #
    ## Initialization
    #
    def _group_by_precursor(self):
        """
        Populate the mapping between precursors and the chromatogram ids.

        The precursor is of type 'PEPT[xx]IDE/3' with an optional (DECOY) tag
        in front. Different modifications will generate different precursors,
        but not different charge states.
        """
        
        openswath_format = self._has_openswath_format(self._run)
        if openswath_format:
            if len( self._run.info['offsets'] ) > 0:
                for key in self._run.info['offsets'].keys():
                    # specific to pymzl, we need to get rid of those two entries
                    if key in ("indexList", "TIC"): continue

                    # The precursor identifier is the second (or third for
                    # decoys) element of the openswath format. Decoys should
                    # get a different identifier!
                    components = key.split("_")
                    trgr_nr = str(components[1])
                    if components[0].startswith("DECOY"):
                        trgr_nr = "DECOY_" + str(components[2])

                    if self._precursor_mapping.has_key(trgr_nr):
                        self._precursor_mapping[trgr_nr].append(key)
                    else:
                        self._precursor_mapping[trgr_nr] = [key]

        else:
            # TODO fallback option!!!
            pass
            raise Exception("Could not parse chromatogram ids ... ")

    def _has_openswath_format(self, run):
        """Checks whether the chromatogram id follows a specific format which
        could allow to map chromatograms to precursors without reading the
        whole file.

        Namely, the format is expected to be [DECOY_]\d*_.*_.* from which one
        can infer that it is openswath format.
        """

        openswath_format = False
        if len( run.info['offsets'] ) > 0:
            keys = run.info['offsets'].keys()
            for key in run.info['offsets'].keys():
                if key in ("indexList", "TIC"): continue
                break

            if len(key.split("_")) >= 3:
                components = key.split("_")
                trgr_nr = components[0]
                if components[0].startswith("DECOY"):
                    trgr_nr = components[1]
                try:
                    trgr_nr = int(trgr_nr)
                    return True
                except ValueError:
                    return False

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
    def get_data_for_transition(self, chrom_id):
        c = self._run[str(chrom_id)] 
        return [ [c.time, c.i] ]

    def get_data_for_precursor(self, precursor):
        """Retrieve data for a specific precursor - data will be as list of
        pairs (timearray, intensityarray)"""

        if not self._precursor_mapping.has_key(str(precursor)):
            return [ [ [0], [0] ] ]

        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            c = self._run[str(chrom_id)] 
            transitions.append([c.time, c.i])

        if len(transitions) == 0: 
            return [ [ [0], [0] ] ]

        return transitions

    def get_range_data(self, precursor):
        return self._range_mapping.get(precursor, [0,0])

    def get_score_data(self, precursor):
        return self._score_mapping.get(precursor, None)

    def get_intensity_data(self, precursor):
        return self._intensity_mapping.get(precursor, None)

    def get_id(self):
        return self._basename

    #
    ## Getters (info)
    #
    def get_transitions_for_precursor(self, precursor):
        return self._precursor_mapping.get(str(precursor), [])

    def get_precursors_for_sequence(self, sequence):
        return self._sequences_mapping.get(sequence, [])

    def get_all_precursor_ids(self):
        return self._precursor_mapping.keys()

    def get_all_peptide_sequences(self):
        return self._sequences_mapping.keys()

class PrecursorModel():

    def __init__(self, chrom_id):
        self.chrom_id = chrom_id

    def getCharge(self):
        try:
            return self.chrom_id.split("/")[1].split("_")[0]
        except Exception:
            return "NA"

    def getFullSequence(self):
        try:
            return self.chrom_id.split("/")[0].split("_")[-1]
        except Exception:
            return "NA"

class DataModel(object):

    def __init__(self):
        self.precursors = set([])
        self.runs = []

    def loadFiles(self, filenames):

        # load new files, clean up ...
        self.runs = []
        self.precursors = set([])

        for f in filenames:
            print "read file", f
            import pymzml
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run = RunDataModel(run_, f)
            self.runs.append(run)
            self.precursors.update(run.get_all_precursor_ids())

    def getStatus(self):
        if len(self.runs) == 0:
            return "Ready"

        tr_cnt = 0
        for r in self.runs:
            tr_cnt += len(r._run.info['offsets']) -2

        return '%s Transitions, %s Peptides (total %s Transitions)' % ( 
            len(self.runs[0]._run.info['offsets']) -2, len(self.runs[0]._sequences_mapping), tr_cnt )

    def get_precursor_list(self):
        return self.precursors

    def get_precursor_tree(self):
        return self._build_tree()

    def _build_tree(self):

        peptide_sequences = set([])
        for r in self.runs:
            peptide_sequences.update( r.get_all_peptide_sequences() )

        ## The peptide sequences are our top-level items
        # print "pepseqs", peptide_sequences
        elements = []
        for seq in peptide_sequences:

            # get all precursors from all runs
            precursors = set([])
            for r in self.runs:
                precursors.update( r.get_precursors_for_sequence(seq) )

            # print "found precursros", precursors
            pelements = []
            for p in precursors:

                # get all transitions from all runs
                transitions = set([])
                for r in self.runs:
                    transitions.update( r.get_transitions_for_precursor(p) )
                tr_elements = []
                pm = PrecursorModel(p)
                for tr in transitions:
                    tr_elements.append(ChromatogramTransition(tr, -1, [], fullName=tr,
                       peptideSequence = pm.getFullSequence(), datatype="Transition") )

                pelements.append(ChromatogramTransition(p, pm.getCharge(), tr_elements, 
                       peptideSequence = pm.getFullSequence(), datatype="Precursor") )
            elements.append(ChromatogramTransition(seq, "NA", pelements, datatype="Peptide", 
                       peptideSequence=pm.getFullSequence()) )
        return elements

    def get_runs(self):
        return self.runs


CHROMTYPES = {
    0 : "Peptide", 
    1 : "Precursor", 
    2 : "Transition"
} 

CHROMTYPES_r = dict([ (v,k) for k,v in CHROMTYPES.iteritems()])

class ChromatogramTransition(object): # your internal structure

    def __init__(self, name, charge, subelements, peptideSequence=None, fullName=None, datatype="Precursor"):
        self.name = name
        self.charge = charge
        self.fullName = fullName
        self.peptideSequence = peptideSequence
        self.subelements = subelements
        self.mytype = CHROMTYPES_r[datatype]

    def getPeptideSequence(self):
        if self.peptideSequence is None:
            return self.name
        return self.peptideSequence

    def getName(self):
        return self.name

    def getType(self):
        return CHROMTYPES[self.mytype]

    def getData(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_data_for_precursor(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_data_for_precursor(prec[0]) 
            else:
                # TODO dont just show the first one!
                pass
                # return run.get_data_for_precursor(prec[0]) 
        elif CHROMTYPES[self.mytype] == "Transition" :
            return run.get_data_for_transition(self.getName()) 
        return [ [ [0], [0] ] ]

    def getRange(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_range_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_range_data(prec[0]) 
        return [ 0,0]

    def getProbScore(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_score_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_score_data(prec[0]) 
        return 1.0

    def getIntensity(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_intensity_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_intensity_data(prec[0]) 
        return -1.0

    def getLabel(self, run):
        # return only the last element
        return [l.split("_")[-1] for l in self._getLabel(run)]

    def _getLabel(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_transitions_for_precursor(self.getName())
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_transitions_for_precursor(prec[0])
            else:
                pass
        elif CHROMTYPES[self.mytype] == "Transition" :
            return [self.getName()]
        return [ "" ]

