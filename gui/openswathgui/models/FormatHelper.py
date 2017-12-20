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

import re

class FormatHelper(object):

    def __init__(self):
        pass

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

    def _has_openswath_format(self, key):
        """
        Checks whether the given key follows a specific format which
        could allow to map chromatograms to precursors without reading the
        whole file.

        Namely, the format is expected to be [DECOY_]\d*_.*_.* from which one
        can infer that it is openswath format.

        possible inputs: 
            DECOY_44736_NVEVIEDDKQGIIR/2_y12
            1002781_TGLC(UniMod:4)QFEDAFTQLSGATPIGAGIDAR_3
        """

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

            charge_is_qualifier = True
            trgr_nr = components[0]
            sequence = components[1]
            if len(sequence.split("/")) > 1:
                sequence = sequence.split("/")[0]
                charge_is_qualifier = False
            qualifier = components[2]

            if not self._is_number(trgr_nr):
                print "Format determination: Could not convert", trgr_nr, "to int."
                return False
            if self._is_number(sequence):
                print "Format determination: Does not look like sequence", sequence
                return False
            elif sequence.find("UniMod") != -1:
                # Looks like sequence, has UniMod tag in it, should be ok
                pass
            else:
                pass
                # valid = re.match('^[A-Z]+$', sequence) is not None
                # if not valid:
                #     print "Format determination: Does not look like sequence", sequence
                #     return False

            if charge_is_qualifier:
                if not self._is_number(qualifier):
                    print "Format determination: Does not look like number", qualifier
                    return False
            elif self._is_number(qualifier):
                print "Format determination: Does look like number", qualifier
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
