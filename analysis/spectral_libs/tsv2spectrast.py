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
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""
from __future__ import print_function

import time, sys
from string import Template
import csv

"""
One sequence of commands to create decoys using spectraST

$ python tsv2spectrast.py library.tsv convertedLibbrary
$ spectrast -cAD -cc -cy1 -c_BIN! -c_ANN -cNconvertedLibrary_plusdecoys convertedLibrary.splib 
$ python sptxt2csv.py convertedLibrary_plusdecoys convertedLibrary_plusdecoys.csv swaths.txt

"""

from msproteomicstoolslib.util import utils
from msproteomicstoolslib.format import speclib_db_lib

infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.reader(open(infile), delimiter='\t')
header = next(reader)
header_d = dict( [ (h,i) for i,h in enumerate(header)] )
lines = list(reader)

assert "group_id" in header
assert "precursor" in header
assert "PeptideSequence" in header
assert "ProteinName" in header
assert "rt_calibrated" in header
assert "charge" in header
assert "product" in header
assert "library_intensity" in header

# group the csv lines by their group id
tr_group = {}
for line in lines:
    if line[header_d["group_id"]] in tr_group:
        tr_group[line[header_d["group_id"]]].append(line)
    else:
        tr_group[line[header_d["group_id"]]] = [line]


# loop through all groups and create a sptxt spectrum for each group
l = speclib_db_lib.Library(1)
cnt = 0
for key in tr_group:
    group = tr_group[key]
    firstline = group[0]
    spectrum = speclib_db_lib.Spectra()
    spectrum.name = firstline[header_d["group_id"]].replace('.', '/')
    spectrum.LibID = 1
    spectrum.MW = firstline[header_d["precursor"]]
    spectrum.precursorMZ = firstline[header_d["precursor"]]
    spectrum.status = "Normal"
    spectrum.full_name = firstline[header_d["group_id"]].replace('.', '/')
    spectrum.sequence = firstline[header_d["PeptideSequence"]]
    spectrum.number_peaks =  len(group)

    ##### The comment

    #AvePrecursorMz=408.0037 BestRawSpectrum=hroest_L120218_005.03628.03628.2 BinaryFileOffset=634 ConsFracAssignedPeaks=0.711 DotConsensus=0.89,0.07;2/48 FracUnassigned=0.25,2/5;0.24,8/20;0.11,13/45 Inst=1/UNKNOWN,48,0 MaxRepSN=327.4 Mods=0 NAA=8 NMC=0 NTT=2 Nreps=48/48 OrigMaxIntensity=1e+05 Parent=407.755 Pep=Tryptic PrecursorIntensity=0 Prob=1.0000 ProbRange=1,1,0.999377,0.9809 Protein=1/H_YeastEnolase ReducedFracIonCurrent=0.361 RepFracAssignedPeaks=0.447 RepNumPeaks=140.6/107.7 RetentionTime=3080.5,3621.3,2704.1 SN=327.4 Sample=1/_mnt_u06_hroest_html_xtandem_michrome_search_xinteract_all_output_xinteract,48,48 Se=1^K48:ex=0.4169/0.8510,fv=6.3620/1.2524,hs=491.9583/70.3180,ns=277.5417/34.1131,pb=0.9994/0.0000 Spec=Consensus TotalIonCurrent=0
    spectrum.comment_dic = {}
    # spectrum.comment_dic["AvePrecursorMz"] = 
    # spectrum.comment_dic["BestRawSpectrum"] = 
    # spectrum.comment_dic["BinaryFileOffset"] = 
    # spectrum.comment_dic["ConsFracAssignedPeaks"] = 
    # spectrum.comment_dic["DotConsensus"] = 
    # spectrum.comment_dic["FracUnassigned"] = 
    # spectrum.comment_dic["Inst"] = 
    # spectrum.comment_dic["MaxRepSN"] = 
    # spectrum.comment_dic["Mods"] = 
    # spectrum.comment_dic["NAA"] = 
    # spectrum.comment_dic["NMC"] = 
    # spectrum.comment_dic["NTT"] = 
    # spectrum.comment_dic["Nreps"] = 
    # spectrum.comment_dic["OrigMaxIntensity"] = 
    # spectrum.comment_dic["Parent"] = 
    # spectrum.comment_dic["Pep"] = 
    # spectrum.comment_dic["PrecursorIntensity"] = 
    # spectrum.comment_dic["Prob"] = 
    # spectrum.comment_dic["ProbRange"] = 
    # spectrum.comment_dic["RepFracAssignedPeaks"] = 
    # spectrum.comment_dic["RepNumPeaks"] = 
    # spectrum.comment_dic["SN"] = 
    # spectrum.comment_dic["Sample"] = 
    # spectrum.comment_dic["Se"] = 
    # spectrum.comment_dic["Spec"] = 
    # spectrum.comment_dic["TotalIonCurrent"] = 
    spectrum.comment_dic["Protein"] = firstline[header_d["ProteinName"]].split("|")[0]
    spectrum.comment_dic["RetentionTime"] = float(firstline[header_d["rt_calibrated"]])
    #spectrum.comment = "Protein=1/%s RetentionTime=%s" % ( firstline[header_d["ProteinName"]].split("|")[0], float(firstline[header_d["rt_calibrated"]])+1000  )
    spectrum.comment = ""
    for k in sorted(spectrum.comment_dic.keys()):
      spectrum.comment += "%s=%s " % (k, spectrum.comment_dic[k])
    spectrum.comment = spectrum.comment[:-1]
    ##### The ptm string
    # | 2|2/7,S,Phospho/-2,K,Methyl| | 
    # 
    # split by "|" and you get: charge, modifications
    # split modifications by "/" and you get:
    #   -- number of modifications
    #   -- modifications themselves
    #
    # then split each mod_string by "," and you get: position, aa, type
    spectrum.ptm_string = "%s|%s|" % (firstline[header_d["charge"]], "0")
    peaks = ""
    for transition in group:
        #peaks += "%s\t%s\t%s\n" % (transition[header_d["product"]], transition[header_d["library_intensity"]], transition[header_d["Annotation"]])
        peaks += "%s\t%s\n" % (transition[header_d["product"]], transition[header_d["library_intensity"]])
        cnt += 1
    spectrum.compress_spectra = peaks
    l.add_spectra(spectrum)

print("wrote ", cnt, "transitions from a total of ", len(tr_group), "peptides")
l.write(outfile)

