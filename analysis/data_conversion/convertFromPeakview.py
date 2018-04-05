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
$Maintainer: Hannes Roest $
$Authors: Hannes Roest $
--------------------------------------------------------------------------
"""

"""
Convert from Peakview format to OpenSWATH format

Example command:
    python convertFromPeakview.py peakview_library.csv openswath_library.csv

"""
import sys,csv

inp = sys.argv[1]
outp = sys.argv[2]

r = csv.reader(open(inp), delimiter="\t")
wr = csv.writer(open(outp, "w"), delimiter="\t")

# write OpenSwath header
wr.writerow([
'PrecursorMz', 'ProductMz', 'Tr_recalibrated', 'transition_name', 'CE', 
'LibraryIntensity', 'transition_group_id', 'ProteinName', 'GroupLabel', 'PeptideSequence', 
'FullUniModPeptideName', 'UniprotID', 'decoy', 'shared', 'confidence', 
'PrecursorCharge', 'FragmentType', 'FragmentCharge', 'FragmentSeriesNumber'])

header = next(r)
# print(header)
hdict = dict( [ (c,i) for i,c in enumerate(header) ] )

already_seen = set([])
tr_cnt = 0
trgr_cnt = 0
if "modification_sequence" in hdict:
    # Ensure we really have a peakview library
    prev_unique_id = ""
    for l in r:

        unique_id = l[hdict["modification_sequence"]] + "/" + l[hdict["prec_z"]]
        if prev_unique_id == "":
            prev_unique_id = unique_id
        if unique_id != prev_unique_id:
            prev_unique_id = unique_id
            trgr_cnt += 1

        already_seen.update( [ l[hdict["modification_sequence"]] ] )
        lnew = [
                l[hdict["Q1"]],
                l[hdict["Q3"]],
                l[hdict["iRT"]],
                str(tr_cnt) + "_" + l[hdict["stripped_sequence"]],
                -1,

                l[hdict["relative_intensity"]],
                str(trgr_cnt) + "_" +  l[hdict["modification_sequence"]] + "/" + l[hdict["prec_z"]],
                # l[hdict["uniprot_id"]],
                l[hdict["UniProtID"]],
                "light",
                l[hdict["stripped_sequence"]],

                l[hdict["modification_sequence"]],
                l[hdict["UniProtID"]],
                0,
                0,
                0,

                l[hdict["prec_z"]],
                l[hdict["frg_type"]],
                l[hdict["frg_z"]],
                l[hdict["frg_nr"]]
        ]
        wr.writerow(lnew)
        tr_cnt += 1
        

