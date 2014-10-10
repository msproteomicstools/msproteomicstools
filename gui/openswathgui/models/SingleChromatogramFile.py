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

class SingleChromatogramFile():
    """Data Model for a single file from one run.

    One run may contain multiple mzML files

    Attributes:
        runid: Current run id

    Private Attributes:
       - _run:        A pymzml.run.Reader object
       - _filename:   Original filename
       - _basename:   Original filename basename
       - _precursor_mapping:   Dictionary { FullPrecursorName : [transition_id, transition_id] }
       - _sequences_mapping:   Dictionary { StrippedSequence : [FullPrecursorName, FullPrecursorName]}

    """

    def __init__(self, run, filename, load_in_memory=False, precursor_mapping = None, sequences_mapping = None):
        import os
        self._run = run
        self._filename = filename
        self._basename = os.path.basename(filename)

        # Map which holds the relationship between a precursor and its
        # corresponding transitions.
        self._precursor_mapping = {}
        # Map which holds the relationship between a sequence and its precursors
        self._sequences_mapping = {}

        self._in_memory = False

        if not precursor_mapping is None and not len(precursor_mapping) == 0 \
          and not sequences_mapping is None and not len(sequences_mapping) == 0:
            self._precursor_mapping = precursor_mapping
            self._sequences_mapping = sequences_mapping
        elif load_in_memory:
            self._group_by_precursor_in_memory()
            self._in_memory = True
            self._group_precursors_by_sequence()
        else:
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
            print "Determined chromatogram identifers to be in openswath formta (11111_PEPTIDE_2), will parse them accordingly."
            
        if openswath_format:
            if len( self._run.info['offsets'] ) > 0:
                for key in self._run.info['offsets'].keys():
                    # specific to pymzl, we need to get rid of those two entries or
                    # when the key is zero
                    if key in ("indexList", "TIC"): continue
                    if len(key) == 0: continue
                    #
                    trgr_nr = self._compute_transitiongroup_from_key(key)
                    tmp = self._precursor_mapping.get(trgr_nr, [])
                    tmp.append(key)
                    self._precursor_mapping[trgr_nr] = tmp

        else:
            print "Fall-back: Could not group chromatograms by their id alone, try using precursor information."
            self._group_by_precursor_by_mass()

    def _group_by_precursor_by_mass(self):
        """
        Populate the mapping between precursors and the chromatogram ids using the chromatogram information.

        Try to use the mass or the precursor tag to infer which chromatograms
        belong to the same peptide. See _compute_transitiongroup_from_precursor
        """

        for key in self._run.info['offsets'].keys():
            # specific to pymzl, we need to get rid of those two entries or
            # when the key is zero
            if key in ("indexList", "TIC"): continue
            if len(key) == 0: continue
            #
            chromatogram = self._run[key]
            if len(chromatogram["precursors"]) == 1:
                precursor = chromatogram["precursors"][0]
                trgr_nr = self._compute_transitiongroup_from_precursor(precursor)
                tmp = self._precursor_mapping.get(trgr_nr, [])
                tmp.append(key)
                self._precursor_mapping[trgr_nr] = tmp

            else:
                print "Something is wrong"
                raise Exception("Could not parse chromatogram ids ... neither ids nor precursors lead to sensible grouping")

    def _group_by_precursor_in_memory(self):
        """
        Same as _group_by_precursor but loads all data in memory.
        """
        openswath_format = self._has_openswath_format(self._run)
        result = {}
        if openswath_format:
            for chromatogram in self._run:
                    key = chromatogram["id"]
                    trgr_nr = self._compute_transitiongroup_from_key(key)
                    tmp = self._precursor_mapping.get(trgr_nr, [])
                    tmp.append(key)
                    self._precursor_mapping[trgr_nr] = tmp

                    import copy
                    cc = copy.copy(chromatogram)
                    result[ key ] = cc
        else:
            raise Exception("Could not parse chromatogram ids ... ")
        
        # we work in memory
        self._run = result

    def _compute_transitiongroup_from_key(self, key):
        """ Transforms an input chromatogram id to a string of [DECOY]_xx/yy

        Possible Input Formats:

        i) [DECOY_]id_xx/yy_zz
        ii) [DECOY_]id_xx_yy

        Where the DECOY_ prefix is optional, xx is the sequence and yy is the
        charge. zz is an additional, optional annotation.

        We want to return only [DECOY]_xx/yy (ensuring that decoys get a
        different identifier than targets).
        """
        components = key.split("_")
        trgr_nr = str(components[1])
        if components[0].startswith("DECOY"):
            trgr_nr = "DECOY_" + str(components[2])

        # Format ii) (second component doesnt contain a slash)
        if trgr_nr.find("/") == -1:
            trgr_nr += "/" + str(components[-1])

        return trgr_nr

    def _compute_transitiongroup_from_precursor(self, precursor):
        """ Transforms an input precursor to a string of [DECOY]_xx/yy
        """
        charge = "0"
        if precursor.has_key("charge"):
            charge = str(precursor["charge"])

        # Check if there is a user-parameter describing the peptide sequence
        # which would allow us to asseble the SEQ/ch pair.
        if precursor.has_key("userParams"):
            if precursor["userParams"].has_key("peptide_sequence"):
                return precursor["userParams"]["peptide_sequence"] + "/" + charge
            if precursor["userParams"].has_key("PeptideSequence"):
                return precursor["userParams"]["PeptideSequence"] + "/" + charge

        if precursor.has_key("mz"):
            # TODO have some reproducible rounding here
            return str(precursor["mz"]) + "GENERIC/" + charge

    def _has_openswath_format(self, run):
        """Checks whether the chromatogram id follows a specific format which
        could allow to map chromatograms to precursors without reading the
        whole file.

        Namely, the format is expected to be [DECOY_]\d*_.*_.* from which one
        can infer that it is openswath format.

        possible inputs: 
            DECOY_44736_NVEVIEDDKQGIIR/2_y12
            1002781_TGLC(UniMod:4)QFEDAFTQLSGATPIGAGIDAR_3
        """

        
        import re
        openswath_format = False
        if len( run.info['offsets'] ) > 0:
            keys = run.info['offsets'].keys()
            for key in run.info['offsets'].keys():
                if key in ("indexList", "TIC"): continue
                if len(key) == 0: continue
                break

            # It is probably not the openswath format with 3 or more entries
            if len(key.split("_")) > 4:
                return False

            if len(key.split("_")) >= 3:
                components = key.split("_")
                # print "split into components", components
                if components[0].startswith("DECOY"):
                    components = components[1:]

                # Exactly 3 components are needed: id, sequence, charge/ion qualifier
                if len(components) != 3:
                    return False

                trgr_nr = components[0]
                sequence = components[1]
                sequence = sequence.split("/")[0]
                qualifier = components[2]

                if not self._is_number(trgr_nr):
                    print "Format determination: Could not convert", trgr_nr, "to int."
                    return False
                if self._is_number(sequence):
                    print "Format determination: Does not look like sequence ", sequence
                    return False
                elif sequence.find("UniMod") != -1:
                    # Looks like sequence, has UniMod tag in it, should be ok
                    pass
                else:
                    valid = re.match('^[A-Z]+$', sequence) is not None
                    if not valid:
                        print "Format determination: Does not look like sequence ", sequence
                        return False

                if self._is_number(qualifier):
                    print "Format determination: Does look like number ", qualifier
                    return False
                
                return True

        # default is false
        return False

    def _is_number(self, n):
        # Sequence should not be a number
        try:
            x = float(n)
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

    def getTransitionCount(self):
        """
        Get total number of transitions 
        """
        if self._in_memory:
            return len(self._run)
        else:
            return len(self._run.info['offsets']) -2

    def get_data_for_precursor(self, precursor):
        """Retrieve raw data for a specific precursor - data will be as list of
        pairs (timearray, intensityarray)"""

        if not self._precursor_mapping.has_key(str(precursor)):
            return [ [ [0], [0] ] ]

        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            chromatogram = self._run[str(chrom_id)] 
            transitions.append([chromatogram.time, chromatogram.i])

        if len(transitions) == 0: 
            return [ [ [0], [0] ] ]

        return transitions

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

    def get_transitions_with_mass_for_precursor(self, precursor):
        """
        Return the transition names prepended with the mass for a specific precursor
        """
        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            chromatogram = self._run[str(chrom_id)] 
            try:
                mz = chromatogram['product']['target_mz']
                if mz > 0.0:
                    transitions.append("%.2f" % mz + " m/z (" + chrom_id + ")")
                else:
                    transitions.append(chrom_id)
            except KeyError:
                transitions.append(chrom_id)
        return transitions

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

