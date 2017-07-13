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
import os, sys
import numpy
import scipy.stats
import random

from msproteomicstoolslib.format.MatrixWriters import getwriter
import msproteomicstoolslib.math.Smoothing as smoothing

if (sys.version_info > (3, 0)):
    xrange = range
else:
    pass

def write_out_matrix_file(matrix_outfile, allruns, multipeptides, fraction_needed_selected,
                          style="none", write_requant=True, aligner_mscore_treshold=1.0):
    matrix_writer = getwriter(matrix_outfile)

    run_ids = [r.get_id() for r in allruns]
    header = ["Peptide", "Protein"]
    for r in allruns:
        fname = "%s_%s" % (os.path.basename(r.orig_filename), r.get_id() )
        header.extend(["Intensity_%s" % fname])
        if style == "RT" or style == 'full':
            header.extend(["RT_%s" % fname])
        if style == "score" or style == 'full':
            header.extend(["score_%s" % fname])

    header.extend(["RT_mean", "RT_std", "pg_pvalue"])

    for i in header:
        matrix_writer.write(i)
    matrix_writer.newline()

    for multipep in sorted(multipeptides, key=lambda x: str(x)):

        # Retrieve all transition group ids available for this precursor group
        # Iterate through all ids and write one line per transition group id
        trgr_ids = set([ trgr.get_id() for prgr in multipep.getPrecursorGroups() for trgr in prgr ])
        for trgr_id in sorted(trgr_ids):

            # Get all selected peakgroups that correspond to the current
            # transition group id and ensure we do not have any twice.
            allruns = [pg.getPeptide().getRunId() for pg in multipep.get_selected_peakgroups() 
                       if pg.getPeptide().get_id() == trgr_id]
            if len(allruns) != len(set(allruns)):
                # TODO test this as well .... 
                raise Exception("Error when writing out matrix, found more than one peakgroup for a run %s" % allruns)

            # Get all selected peakgroups that correspond to the current transition group id 
            selected_peakgroups = dict([(pg.getPeptide().getRunId(), pg) 
                for pg in multipep.get_selected_peakgroups() if pg.getPeptide().get_id() == trgr_id])

            # Skip empty lines or lines that have too few entries
            if len(selected_peakgroups) == 0:
                continue
            if (len(selected_peakgroups) * 1.0 / len(allruns) < fraction_needed_selected): 
                continue

            # Write first two columns of the matrix
            for i in [trgr_id, multipep.find_best_peptide_pg().getPeptide().getProteinName()]:
                matrix_writer.write(i)

            # Write other columns (one or two per run, depending on format)
            rts = []
            for rid in run_ids:
                pg = selected_peakgroups.get(rid, None)

                if not write_requant:
                    if not pg is None and pg.get_fdr_score() > 1.0:
                        pg = None

                if pg is None:
                    matrix_writer.write("NA")
                    if style == "RT" or style == "full":
                        matrix_writer.write("NA")
                    if style == "score" or style == "full":
                        matrix_writer.write("NA")
                else:
                    if pg.get_fdr_score() > 1.0:
                        color = 'r'
                    elif pg.get_fdr_score() > aligner_mscore_treshold:
                        color = 'b'
                    else:
                        color = 'd'

                    matrix_writer.write(pg.get_intensity(), color=color)

                    if style == "RT" or style == "full":
                        matrix_writer.write(pg.get_normalized_retentiontime(), color=color)
                    if style == "score" or style == "full":
                        matrix_writer.write(pg.get_fdr_score(), color=color)

                if not pg is None:
                    rts.append(pg.get_normalized_retentiontime())

            # The d_score is a z-score which computed on the null / decoy
            # distribution which is (assumed) gaussian with u = 0, sigma = 1
            # -> we thus compute a p-value from the z-score and assuming
            # independent measurements, we multiply the p-values to compute a
            # peakgroup p-value.
            # We use norm.sf (1-cdf) on the vector of z-scores.
            pvals = [float(pg.get_dscore()) for k,pg in selected_peakgroups.items() if
                     not pg is None and not pg.get_dscore() is None]
            pvalue = numpy.prod(scipy.stats.norm.sf(pvals))

            for i in [numpy.mean(rts), numpy.std(rts), pvalue]:
                matrix_writer.write(i)
            matrix_writer.newline()

    del matrix_writer

def addDataToTrafo(tr_data, run_0, run_1, spl_aligner, multipeptides,
                   realign_method, max_rt_diff, topN=5, sd_max_data_length=5000, force=False):
    id_0 = run_0.get_id()
    id_1 = run_1.get_id()

    if id_0 == id_1:
        null = smoothing.SmoothingNull()
        tr_data.addTrafo(id_0, id_1, null)
        tr_data.addTrafo(id_1, id_0, null)
        return

    # Data
    data_0, data_1 = spl_aligner._getRTData(run_0, run_1, multipeptides)
    tr_data.addData(id_0, data_0, id_1, data_1)

    # import pylab
    # pylab.scatter(data_0, data_1)
    # pylab.savefig('data_%s_%s.pdf' % (run_0, run_1) )
    # pylab.clf()
    # pylab.scatter(data_0, data_1)
    # pylab.xlim(2300, 2600)
    # pylab.ylim(2300, 2600)
    # pylab.savefig('data_%s_%s_zoom.pdf' % (run_0, run_1) )
    # pylab.clf()

    if len(data_0) == 0:
        print("Warning, zero data!")
        if force:
            null = smoothing.SmoothingNull()
            tr_data.addTrafo(id_0, id_1, null)
            tr_data.addTrafo(id_1, id_0, null)
            return
        else:
            raise Exception("No data available for alignment %s vs %s" % (id_0, id_1) )

    # Smoothers
    sm_0_1 = smoothing.getSmoothingObj(realign_method, topN=topN,
                                       max_rt_diff=max_rt_diff,
                                       min_rt_diff=0.1, removeOutliers=False,
                                       tmpdir=None)
    sm_1_0 = smoothing.getSmoothingObj(realign_method, topN=topN,
                                       max_rt_diff=max_rt_diff,
                                       min_rt_diff=0.1, removeOutliers=False,
                                       tmpdir=None)
    # Initialize smoother
    sm_0_1.initialize(data_0, data_1)
    sm_1_0.initialize(data_1, data_0)

    # Compute error for alignment (standard deviation)
    stdev_0_1 = 0.0
    stdev_1_0 = 0.0
    if sd_max_data_length > 0:
        sample_idx = random.sample( xrange(len(data_0)), min(sd_max_data_length, len(data_0))  )
        data_0_s = [data_0[i] for i in sample_idx]
        data_1_s = [data_1[i] for i in sample_idx]
        data0_aligned = sm_0_1.predict(data_0_s)
        stdev_0_1 = numpy.std(numpy.array(data_1_s) - numpy.array(data0_aligned))
        data1_aligned = sm_1_0.predict(data_1_s)
        stdev_1_0 = numpy.std(numpy.array(data_0_s) - numpy.array(data1_aligned))
        print("stdev for", id_0, id_1, stdev_0_1, " / ", stdev_1_0, "on data length", len(data_0_s))

    # Add data and trafo description.
    # The CyLightTransformationData actually requires to get a specific type of
    # transformation, the CyLinearInterpolateWrapper which may not be directly
    # passed to this function. We will try to recover the underlying linear
    # wrapper and then stick it into the tr_data object. If this fails, we just
    # revert to the regular behavior.
    try:
        sm_0_1_lwp = sm_0_1.internal_interpolation.getLWP()
        sm_1_0_lwp = sm_1_0.internal_interpolation.getLWP()
        tr_data.addTrafo(id_0, id_1, sm_0_1_lwp, stdev_0_1)
        tr_data.addTrafo(id_1, id_0, sm_1_0_lwp, stdev_1_0)
    except Exception:
        tr_data.addTrafo(id_0, id_1, sm_0_1, stdev_0_1)
        tr_data.addTrafo(id_1, id_0, sm_1_0, stdev_1_0)

