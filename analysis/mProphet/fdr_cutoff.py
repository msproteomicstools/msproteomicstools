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
filename = sys.argv[1]
filename_out = sys.argv[2]
sort_by = sys.argv[3]
decoy_name = sys.argv[4]
fdr_cutoff = float(sys.argv[5])
reverse = sys.argv[6]
reverse = (reverse == 'TRUE')
use_fdr_score = sys.argv[7]
use_fdr_score = (use_fdr_score == 'TRUE')

reader = csv.reader(open(filename), delimiter='\t')
header = next(reader)
header_d = dict( [ (h,i) for i,h in enumerate(header)] )
try:
  sort_pos = header_d[sort_by]
except KeyError:
    print("Could not find key", sort_by, "in", header_d)
    sys.exit()


print("Sort by '%s' (position %s) and reverse %s" % (sort_by, sort_pos, reverse))
lines = list(reader)
lines.sort(key=lambda x: float(x[sort_pos]), reverse=reverse)

# try to estimate an empirical FDR and only keep those lines above the fdr
decoys_ = 0
for i,line in enumerate(lines):
    is_decoy = False
    for cell in line: 
        if cell.find(decoy_name) != -1: 
            is_decoy = True

    if is_decoy:
        decoys_ += 1

    if use_fdr_score:
      fdr = float(line[sort_pos])
    elif (i+1-decoys_) > 0:
      fdr = decoys_ * 1.0/ (i+1-decoys_)
    else: fdr = 1.0
    line.append(fdr)
    if(fdr > fdr_cutoff): break


print("Stop after %s entries of which %s are decoys (est. fdr=%s)" % (i, decoys_, fdr))
above_cutoff = lines[:i]

header.append('FDREstimate')
writer = csv.writer(open(filename_out, 'w'), delimiter='\t')
writer.writerow(header)
writer.writerows(above_cutoff)
