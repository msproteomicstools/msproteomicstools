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
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""

import sys

from msproteomicstoolslib.format import pepXMLReader

infile = sys.argv[1]
outfile = sys.argv[2]
r = pepXMLReader.pepXMLReader(infile)
import csv
writer = csv.writer(open(outfile, 'w'))
writer.writerow(
    ["Spectrum", "Scan number", "Peptide_Sequence", "PrecursorMass", "Massdiff", "Expect value", "Pvalue", "MatchedIons", "TotalIons",  "Protein"])
for hit in r.parse_all():
  #exp = float(hit.expect)
  #pval = float(hit.pvalue)
  #exp = str(hit.expect).replace('e', 'E')
  pval = "NA"
  if hasattr(hit, "pvalue"): pval = "%.20f"  % hit.pvalue
  exp = "NA"
  if hasattr(hit, "expect"): exp = "%.20f"  % hit.expect
  writer.writerow(
    [hit.spectrum_query.spectrum, hit.scan_number, hit.peptide,
     hit.spectrum_query.precursorMass, hit.massdiff, exp, pval,
     hit.matched_ions, hit.total_ions, hit.protein_descr]
  )

