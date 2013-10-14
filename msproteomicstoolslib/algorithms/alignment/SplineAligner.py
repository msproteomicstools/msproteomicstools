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

import numpy
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection

class SplineAligner():
    """
    Use the datasmoothing part of msproteomicstoolslib to align 2 runs in
    retention times using splines.


    >>> spl_aligner = SplineAligner()
    >>> transformations = spl_aligner.rt_align_all_runs(this_exp, multipeptides, options.alignment_score, options.use_scikit)
    """
    def __init__(self):
      self.transformation_collection = TransformationCollection()

    def _determine_best_run(self, experiment, alignment_fdr_threshold):

        maxcount = -1
        bestrun = -1
        for run in experiment.runs:
            cnt = 0
            for peptide in run.all_peptides.values():
                if peptide.get_decoy(): continue
                pg = peptide.get_best_peakgroup()
                if pg.get_fdr_score() < alignment_fdr_threshold:
                    cnt += 1
            if cnt > maxcount:
                maxcount = cnt
                bestrun = run.get_id()
        print "Found best run", bestrun, "with %s features above the cutoff of %s%%" % (maxcount, alignment_fdr_threshold)
        return [r for r in experiment.runs if r.get_id() == bestrun][0]

    def _spline_align_runs(self, bestrun, run, multipeptides, alignment_fdr_threshold, use_scikit, use_linear):
        """Will align run against bestrun"""

        # get those peptides we want to use for alignment => for this use the mapping
        # data1 = reference data (master)
        # data2 = data to be aligned (slave)
        data1 = []
        data2 = []
        for m in multipeptides:
            try: 
                ref_pep = m.get_peptide(bestrun.get_id()).get_best_peakgroup()
                align_pep = m.get_peptide(run.get_id()).get_best_peakgroup()
            except KeyError: 
                # it is possible that for some, no peak group exists in this run
                continue
            if ref_pep.peptide.get_decoy() or align_pep.peptide.get_decoy(): continue
            if ref_pep.get_fdr_score() < alignment_fdr_threshold and align_pep.get_fdr_score() < alignment_fdr_threshold:
                data1.append(ref_pep.get_normalized_retentiontime())
                data2.append(align_pep.get_normalized_retentiontime())

        # from run to bestrun
        self.transformation_collection.addTransformationData([data2, data1], run.get_id(), bestrun.get_id() )

        all_pg = []
        for pep in run.all_peptides.values():
            all_pg.extend( [ (pg.get_normalized_retentiontime(), pg.get_feature_id()) for pg in pep.get_all_peakgroups()] )

        rt_eval = [ pg[0] for pg in all_pg]

        # Since we want to predict how to convert from slave to master, slave
        # is first and master is second.
        sm = smoothing.get_smooting_operator(use_scikit = use_scikit, use_linear = use_linear)
        sm.initialize(data2, data1)
        aligned_result = sm.predict(rt_eval)

        data2_aligned = sm.predict(data2)
        self.transformation_collection.addTransformedData(data2_aligned, run.get_id(), bestrun.get_id() )

        # The two methods produce very, very similar results
        # but R is faster => prefer to use R when possible.
        # hist(aligned_result - aligned_result_2, 100)
        # numpy.std(aligned_result - aligned_result_2)
        # 0.66102016517870454
        # numpy.median(aligned_result - aligned_result_2)
        # -0.020456989235640322

        print "Will align run %s against %s, using %s features" % (run.get_id(), bestrun.get_id(), len(data1))
        print "  Computed stdev", numpy.std(numpy.array(data1) - numpy.array(data2_aligned)), \
                  "and median", numpy.median(numpy.array(data1) - numpy.array(data2_aligned))

        # now re-populate the peptide data!
        i = 0
        for pep in run.all_peptides.values():
            mutable = [list(pg) for pg in pep.peakgroups_]
            for k in range(len(mutable)):
                mutable[k][2] = aligned_result[i]
                i += 1
            pep.peakgroups_ = [ tuple(m) for m in mutable]

    def rt_align_all_runs(self, experiment, multipeptides, alignment_fdr_threshold = 0.0001, use_scikit=False, use_linear=False):
        """ Align all runs contained in an MRExperiment

        Args:
            experiment(MRExperiment): a collection of runs
            multipeptides(list(multipeptides)): a list of Multipeptide derived from the above expriment
        """

        print "Will re-align runs"
        # spl_aligner = SplineAligner()

        # get the best run (e.g. the one with the most ids below threshold)
        bestrun = self._determine_best_run(experiment, alignment_fdr_threshold)

        ## spl_aligner.transformation_collection = experiment.transformation_collection
        self.transformation_collection.setReferenceRunID( bestrun.get_id() )

        # go through all runs and align two runs at a time
        for run in experiment.runs:
            if run.get_id() == bestrun.get_id(): continue # do not align reference run itself
            self._spline_align_runs(bestrun, run, multipeptides, alignment_fdr_threshold, use_scikit, use_linear)

        return self.transformation_collection

