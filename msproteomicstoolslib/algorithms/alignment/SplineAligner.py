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
import numpy
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection

class TransformationError():
    """
    Class to store the error of a given transformation
    """
    def __init__(self):
        self.transformations = {}

    def getStdev(self):
        """
        Retrieve the standard deviation as a measurement of the error of a transformation
        """
        for tr1 in self.transformations.values():
            for tr2 in tr1.values():
                yield tr2[0]

class SplineAligner():
    """
    Use the datasmoothing part of msproteomicstoolslib to align two runs in
    retention times using splines.

    >>> spl_aligner = SplineAligner()
    >>> transformations = spl_aligner.rt_align_all_runs(this_exp, multipeptides, options.alignment_score, options.use_scikit)
    """
    def __init__(self, alignment_fdr_threshold = 0.0001, smoother="lowess", external_r_tmpdir=None, maxdata=-1, experiment=None):
      self.transformation_collection = TransformationCollection()
      self.alignment_fdr_threshold_ = alignment_fdr_threshold
      self.smoother = smoother
      self.tmpdir_ = external_r_tmpdir
      self.max_data_ = maxdata
      self._cacher = None
      self._cy_cacher = None
      self._experiment = experiment

      try:
          from msproteomicstoolslib.algorithms.alignment.DataCacher import CyDataCacher
          self._cy_cacher = CyDataCacher()
      except ImportError:
          print("WARNING: cannot import CyDataCacher, will use Python version (slower).")

    def _determine_best_run(self, experiment):

        maxcount = -1
        bestrun = -1
        for run in experiment.runs:
            cnt = 0
            for prgroup in run:
                for peptide in prgroup:
                    if peptide.get_decoy(): continue
                    pg = peptide.get_best_peakgroup()
                    if pg.get_fdr_score() < self.alignment_fdr_threshold_:
                        cnt += 1
            if cnt > maxcount:
                maxcount = cnt
                bestrun = run.get_id()
        print("Found best run", bestrun, "with %s features above the cutoff of %s%%" % (maxcount, self.alignment_fdr_threshold_))
        return [r for r in experiment.runs if r.get_id() == bestrun][0]

    def _getRTData(self, bestrun, run, multipeptides):
        """ Return retention time data for reference and slave run """

        if self._experiment is not None:
            return self._getRTData_cached(bestrun, run, multipeptides)
        else:
            return self._getRTData_legacy(bestrun, run, multipeptides)

    def _cache_RT_data(self, bestrun, run, multipeptides):

        self._cacher = []
        for m in multipeptides:
            cached_vals = []
            is_decoy = False
            for r in self._experiment.runs:
                val = None
                if m.hasPrecursorGroup(r.get_id()):

                    al_pg = [pg for pg in m.getPrecursorGroup(r.get_id()).getAllPeakgroups()
                                   if pg.get_fdr_score() < self.alignment_fdr_threshold_]

                    # We need to have a single, good peak group below the threshold (not a decoy)
                    if len(al_pg) == 1:
                        pep = m.getPrecursorGroup(r.get_id()).getOverallBestPeakgroup()
                        if not pep.getPeptide().get_decoy() and pep.get_fdr_score() < self.alignment_fdr_threshold_:
                            val = (pep.get_fdr_score(), pep.get_normalized_retentiontime())

                cached_vals.append(val)

            # only append with at least 2 values ...
            if len([v for v in cached_vals if not v is None] ) > 1:
                if self._cy_cacher is not None:
                    self._cy_cacher.appendValuesForPeptide(cached_vals)
                else:
                    self._cacher.append(cached_vals)

    def _getRTData_cached(self, bestrun, run, multipeptides):
        """ Return retention time data for reference and slave run """

        if self._cacher is None:
            self._cache_RT_data(bestrun, run, multipeptides)

        run_nr = [k for k,r in enumerate(self._experiment.runs) if r.get_id() == run.get_id() ][0]
        bestrun_nr = [k for k,r in enumerate(self._experiment.runs) if r.get_id() == bestrun.get_id() ][0]

        if self._cy_cacher is not None:

            return self._cy_cacher.retrieveValues(bestrun_nr, run_nr)

        data_tmp = []
        for m in self._cacher:

            bestrund = m[ bestrun_nr ] # data1
            rund = m[ run_nr ] # data2

            # Skip empty entries
            if rund is None or bestrund is None:
                continue

            data_tmp.append( (
                min( rund[0], bestrund[0]),
                bestrund[1], rund[1]) ) # data1, data2

        maxdata = self.max_data_
        if maxdata == -1:
            # -1 means take all data
            maxdata = len(data_tmp)

        data1 = []
        data2 = []
        for fdr, d1, d2 in sorted(data_tmp)[:maxdata]:
            data1.append(d1)
            data2.append(d2)

        return data1, data2

    def _getRTData_legacy(self, bestrun, run, multipeptides):
        """ Return retention time data for reference and slave run """

        # data1 = reference data (master)
        # data2 = data to be aligned (slave)
        data1 = []
        data2 = []

        data_tmp = []
        cnt_multiple = 0

        for m in multipeptides:

            try:
                len_ali = len([pg for pg in m.getPrecursorGroup(run.get_id()).getAllPeakgroups()
                               if pg.get_fdr_score() < self.alignment_fdr_threshold_])
                len_ref = len([pg for pg in m.getPrecursorGroup(bestrun.get_id()).getAllPeakgroups()
                               if pg.get_fdr_score() < self.alignment_fdr_threshold_])

                # Do not consider peakgroups that are missing in one run
                # Do not consider peakgroups that have more than one good peakgroup
                if len_ali != 1 or len_ref != 1:
                    if len_ali > 1 or len_ref > 1:
                        cnt_multiple += 1
                    continue

                ref_pep = m.getPrecursorGroup(bestrun.get_id()).getOverallBestPeakgroup()
                align_pep = m.getPrecursorGroup(run.get_id()).getOverallBestPeakgroup()
            except KeyError:
                # it is possible that for some, no peak group exists in this run
                continue

            # Do not use decoy peptides
            if ref_pep.peptide.get_decoy() or align_pep.peptide.get_decoy():
                continue

            if ref_pep.get_fdr_score() < self.alignment_fdr_threshold_ and \
               align_pep.get_fdr_score() < self.alignment_fdr_threshold_:

                # data1.append(ref_pep.get_normalized_retentiontime())
                # data2.append(align_pep.get_normalized_retentiontime())
                data_tmp.append( (
                    ref_pep.get_fdr_score(),
                    ref_pep.get_normalized_retentiontime(),
                    align_pep.get_normalized_retentiontime()
                ) )

        if cnt_multiple > len(multipeptides) * 0.8 :
            print ("")
            print ("  Warning: Most of your data has more than one peakgroup with a score better than %s."  % self.alignment_fdr_threshold_)
            print ("  This may be a problem for the alignment, please consider adjusting the --alignment_score option."  )

        maxdata = self.max_data_
        if maxdata == -1:
            # -1 means take all data
            maxdata = len(data_tmp)

        for fdr, d1, d2 in sorted(data_tmp)[:maxdata]:
            data1.append(d1)
            data2.append(d2)

        return data1,data2

    def _spline_align_runs(self, bestrun, run, multipeptides):
        """Will align run against bestrun"""

        sm = smoothing.getSmoothingObj(smoother = self.smoother, tmpdir = self.tmpdir_)

        # get those peptides we want to use for alignment => for this use the mapping
        # data1 = reference data (master)
        # data2 = data to be aligned (slave)
        data1,data2 = self._getRTData(bestrun, run, multipeptides)

        if len(data2) < 2:
            print("No common identifications between %s and %s. Only found %s features below a cutoff of %s" % ( 
                run.get_id(), bestrun.get_id(), len(data1), self.alignment_fdr_threshold_) )
            print("If you ran the feature_alignment.py script, try to skip the re-alignment step (e.g. remove the --realign_runs option)." )
            raise Exception("Not enough datapoints (less than 2 datapoints).")

        # Since we want to predict how to convert from slave to master, slave
        # is first and master is second.
        sm.initialize(data2, data1)
        data2_aligned = sm.predict(data2)

        # Store transformation in collection (from run to bestrun)
        self.transformation_collection.addTransformationData([data2, data1], run.get_id(), bestrun.get_id() )
        self.transformation_collection.addTransformedData(data2_aligned, run.get_id(), bestrun.get_id() )

        stdev = numpy.std(numpy.array(data1) - numpy.array(data2_aligned))
        median = numpy.median(numpy.array(data1) - numpy.array(data2_aligned))
        print("Will align run %s against %s, using %s features" % (run.get_id(), bestrun.get_id(), len(data1)) )
        print("  Computed stdev", stdev, "and median", median )

        # Store error for later
        d = self.transformation_error.transformations.get(run.get_id(), {})
        d[bestrun.get_id()] = [stdev, median]
        self.transformation_error.transformations[ run.get_id() ] = d

        # Now predict on _all_ data and write this back to the data
        i = 0
        all_pg = []
        for prgr in run:
            for pep in prgr:
                all_pg.extend( [ (pg.get_normalized_retentiontime(), pg.get_feature_id()) for pg in pep.get_all_peakgroups()] )
        rt_eval = [ pg[0] for pg in all_pg]
        aligned_result = sm.predict(rt_eval)
        for prgr in run:
            for pep in prgr:
                # TODO hack -> direct access to the internal peakgroups object
                mutable = [list(pg) for pg in pep.peakgroups_]
                for k in range(len(mutable)):
                    mutable[k][2] = aligned_result[i]
                    i += 1
                pep.peakgroups_ = [ tuple(m) for m in mutable]

    def rt_align_all_runs(self, experiment, multipeptides):
        """ Align all runs contained in an MRExperiment

        Args:
            experiment(MRExperiment): a collection of runs
            multipeptides(list(multipeptides)): a list of Multipeptide derived from the above expriment
        """

        print("Will re-align runs" )

        # get the best run (e.g. the one with the most ids below threshold)
        bestrun = self._determine_best_run(experiment)

        ## spl_aligner.transformation_collection = experiment.transformation_collection
        self.transformation_collection.setReferenceRunID( bestrun.get_id() )
        self.transformation_error = TransformationError()

        # go through all runs and align two runs at a time
        for run in experiment.runs:
            if run.get_id() == bestrun.get_id(): continue # do not align reference run itself
            self._spline_align_runs(bestrun, run, multipeptides)

        return self.transformation_collection

    def getTransformationError(self):
        """
        Get the error of the transformation

        Returns:
            transformation_error(:class:`.TransformationError`) : the error of the transformation
        """
        return self.transformation_error

