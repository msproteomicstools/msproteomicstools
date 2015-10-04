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

import csv, sys

# This program counts how many precursors, peptides and proteins are present in
# a csv table. This can be useful after an mProphet run to simply look at the
# result. Please make sure that decoys are removed!

filename = sys.argv[1]

reader = csv.reader(open(filename), delimiter='\t')
header = next(reader)
header_d = dict( [ (h,i) for i,h in enumerate(header)] )

#print header
lines = list(reader)

# try to estimate an empirical FDR and only keep those lines above the fdr
decoys_ = 0
peptides = {}
proteins = {}
pep_per_prot = {}
for i,line in enumerate(lines):
    peptides[ line[header_d['Sequence'] ]] = 0
    proteins[ line[header_d['ProteinName'] ]] = 0
    if line[header_d['ProteinName'] ] not in pep_per_prot:
      pep_per_prot[  line[header_d['ProteinName'] ] ] = {}
    pep_per_prot[  line[header_d['ProteinName'] ] ][ line[header_d['Sequence'] ]] = 0

single_hits = len([0 for k,v in pep_per_prot.items() if len(v) == 1])
print("Found precursors" , i, " and peptides", len(peptides), " and proteins", len(proteins), "(single hits", single_hits, "multiple hits", len(proteins) - single_hits, ")")




