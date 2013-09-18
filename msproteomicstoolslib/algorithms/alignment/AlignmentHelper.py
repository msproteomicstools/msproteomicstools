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

import csv, os
import numpy

def write_out_matrix_file(matrix_outfile, allruns, multipeptides, fraction_needed_selected):
    import scipy.stats
    matrix_writer = csv.writer(open(matrix_outfile, "w"), delimiter="\t")
    run_ids = [r.get_id() for r in allruns]
    header = ["Peptide", "Protein"]
    for r in allruns:
        fname = "%s_%s" % (os.path.basename(r.orig_filename), r.get_id() )
        header.extend(["Intensity_%s" % fname, "RT_%s" % fname])
    header.extend(["RT_mean", "RT_std", "pg_pvalue"])
    matrix_writer.writerow(header)
    for m in multipeptides:
        line = [m.get_id(), m.find_best_peptide_pg().peptide.protein_name]
        rts = []
        selected_peakgroups = m.get_selected_peakgroups()
        if (len(selected_peakgroups)*1.0 / len(allruns) < fraction_needed_selected) : continue
        for rid in run_ids:
            pg = None
            if m.has_peptide(rid):
                pg = m.get_peptide(rid).get_selected_peakgroup()
            if pg is None:
                line.extend(["NA", "NA"])
            else:
                line.extend([pg.get_intensity(), pg.get_normalized_retentiontime()])
                rts.append(pg.get_normalized_retentiontime())

        # The d_score is a z-score which computed on the null / decoy
        # distribution which is (assumed) gaussian with u = 0, sigma = 1 
        # -> we thus compute a p-value from the z-score and assuming
        # independent measurements, we multiply the p-values to compute a
        # peakgroup p-value.
        # We use norm.sf (1-cdf) on the vector of z-scores.
        pvals = [ pg.get_dscore() for pg in m.get_selected_peakgroups() if not pg is None and not pg.get_dscore() is None]
        pvalue = numpy.prod(scipy.stats.norm.sf(pvals))
        line.extend([numpy.mean(rts), numpy.std(rts), pvalue ])
        matrix_writer.writerow(line)
    del matrix_writer

