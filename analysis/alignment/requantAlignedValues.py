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
from __future__ import division

import os, sys, csv, time
import numpy
import argparse
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.data_structures.Precursor import GeneralPrecursor, Precursor
from msproteomicstoolslib.data_structures.PrecursorGroup import PrecursorGroup
import msproteomicstoolslib.format.TransformationCollection as transformations
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
from msproteomicstoolslib.algorithms.alignment.AlignmentMST import getDistanceMatrix
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree
import msproteomicstoolslib.algorithms.alignment.AlignmentHelper as helper
from msproteomicstoolslib.algorithms.alignment.BorderIntegration import \
        integrationBorderShortestPath, integrationBorderShortestDistance, integrationBorderReference
import msproteomicstoolslib.math.Smoothing as smoothing
from feature_alignment import Experiment

# The Window overlap which needs to be taken into account when calculating from which swath window to extract!
SWATH_EDGE_SHIFT = 1
VERBOSE = False

class ImputeValuesHelper(object):
    """
    Static object with some helper methods.
    """

    @staticmethod
    def select_correct_swath(swath_chromatograms, mz):
        """Select the correct chromatogram

        Args:
            swath_chromatograms(dict): containing the objects pointing to the original chrom mzML (see runImputeValues)
            mz(float): the mz value of the precursor
        """
        mz = mz + SWATH_EDGE_SHIFT
        swath_window_low = int(mz // 25) * 25
        swath_window_high = int(mz // 25) * 25 + 25
        res = {}
        for k, v in swath_chromatograms.items():
            # TODO smarter selection here
            selected = [vv for prec_mz, vv in v.items() if prec_mz >= swath_window_low and prec_mz < swath_window_high]
            if len(v) == 1: 
                # We have merged chrom.mzML file (only one file)
                selected = list(v.values())
            if len(selected) == 1: 
                res[k] = selected[0]
        return res

class SwathChromatogramRun(object):
    """ A single SWATH LC-MS/MS run.

    Each run may contain multiple files (split up by swath).
    """

    def __init__(self):
        self.chromfiles = []

    def parse(self, runid, files):
        """ Parse a set of files which all belong to the same experiment 
        """
        for i,f in enumerate(files):
            import pymzml
            if f.find("rtnorm") != -1 or f.find("ms1scan") != -1: 
                continue
            run = pymzml.run.Reader(f, build_index_from_scratch=True)
            if not run.info["seekable"]:
                print("Could not open mzML file %s in seekable mode. Please use unzipped files, preferably with an index (using the indexedmzML standard).") 
                raise Exception("Could not properly read file", f)
            self.chromfiles.append(run)

    def getChromatogram(self, chromid):
        for cfile in self.chromfiles:
            if chromid in cfile.info["offsets"]:
                return cfile[chromid]
        return None

class SwathChromatogramCollection(object):
    """ A collection of multiple SWATH LC-MS/MS runs.

    Each single run is represented as a SwathChromatogramRun and accessible
    through a run id.

    >>> mzml_files = ["file1.mzML", "file2.mzML"]
    >>> runMapping = {"file1.mzML": "run1", "file2.mzML" : "run2"}
    >>> chromatograms = SwathChromatogramCollection()
    >>> chromatograms.parseFromMzML(mzml_files, runMapping)

    >>> chromatogram = chromatograms.getChromatogram("run1", "ChromatogramId")
    >>> chromatogram = chromatograms.getChromatogram("run2", "ChromatogramId")
    """

    def __init__(self):
        self.allruns = {}
        self.cached_run = None
        self.cache = {}

    def createRunCache(self, runid):
        self.cache = {}
        self.cached_run = runid

        import copy
        for run in self.allruns[runid].chromfiles:
            for chromid, value in run.info['offsets'].items():
              if value is None: continue
              self.cache[ chromid ] = copy.copy(run[chromid])

    def _getChromatogramCached(self, runid, chromid):
        assert runid == self.cached_run
        return self.cache.get( chromid, None)

    def getChromatogram(self, runid, chromid):
        if runid == self.cached_run:
            return self._getChromatogramCached(runid, chromid)

        if not runid in self.allruns:
            return None
        return self.allruns[runid].getChromatogram(chromid)

    def parseFromTrafoFiles(self, trafo_fnames):
        """ Parse a set of different experiments from the .tr files

        The mzML files belonging to the same run are assumed to be in the same
        folder as the .tr files.
        """
        swath_chromatograms = {}
        for filename in trafo_fnames:
            # get the run id
            start = time.time()
            f = open(filename, "r")
            header = next(f).split("\t")
            f.close()
            all_swathes = {}
            runid = header[1]
            import glob
            dname = os.path.dirname(filename)
            files = glob.glob(os.path.join(dname + "/*.mzML") )

            swathrun = SwathChromatogramRun()
            swathrun.parse(runid, files)
            self.allruns[runid] = swathrun
            print("Parsing chromatograms in", filename, "took %0.4fs" % (time.time() - start))

    def parseFromMzML(self, mzML_files, runIdMapping):
        """ Parse a set of different experiments.

        Args:
            mzML_files(list(filename)): a list of filenames of the mzML files
            runIdMapping(dict): a dictionary mapping each filename to a run id
        """
        swath_chromatograms = {}
        for filename in mzML_files:
            start = time.time()
            runid = runIdMapping[filename]
            swathrun = SwathChromatogramRun()
            swathrun.parse(runid, [filename])
            self.allruns[runid] = swathrun
            print("Parsing chromatograms in", filename, "took %0.4fs" % (time.time() - start))

    def getRunIDs(self):
        return list(self.allruns.keys())

# Main entry points:
# runSingleFileImputation -> for tree based, new imputation
# runImputeValues -> old method using .tr files

def runSingleFileImputation(options, peakgroups_file, mzML_file, method, is_test):
    """Impute values across chromatograms

    Args:
        peakgroups_file(filename): CSV file containing all peakgroups
        mzML_file(filename): mzML file containing chromatograms
        method(string): which method to use for imputation ("singleShortestPath", "singleClosestRun")
        is_test(bool): whether test mode should be used
    Returns:
        A tuple of:
            new_exp(AlignmentExperiment): experiment containing the aligned peakgroups
            multipeptides(list(AlignmentHelper.Multipeptide)): list of multipeptides
            rid(string): run id of the analyzed run

    This function will read the csv file with all peakgroups as well as the
    provided chromatogram file (.chrom.mzML). It will then try to impute
    missing values for those peakgroups where no values is currently present,
    reading the raw chromatograms.
    """

    # We do not want to exclude any peakgroups for noiseIntegration (we assume
    # that alignment has already happened)
    fdr_cutoff_all_pg = 1.0

    start = time.time()
    reader = SWATHScoringReader.newReader([peakgroups_file],
                                          options.file_format,
                                          readmethod="complete",
                                          enable_isotopic_grouping = not options.disable_isotopic_grouping,
                                          read_cluster_id=True)
    new_exp = Experiment()
    new_exp.runs = reader.parse_files()
    multipeptides = new_exp.get_all_multipeptides(fdr_cutoff_all_pg, verbose=False)
    print("Parsing the peakgroups file took %ss" % (time.time() - start) )

    mapping = {}
    precursors_mapping = {}
    sequences_mapping = {}
    protein_mapping = {}
    inferMapping([ mzML_file ], [ peakgroups_file ], mapping, precursors_mapping, sequences_mapping, protein_mapping, verbose=False)
    mapping_inv = dict((v[0], k) for k, v in mapping.items())
    if VERBOSE:
        print (mapping)

    # Do only a single run : read only one single file
    start = time.time()
    swath_chromatograms = SwathChromatogramCollection()
    swath_chromatograms.parseFromMzML([ mzML_file ], mapping_inv)
    print("Reading the chromatogram files took %ss" % (time.time() - start) )
    assert len(swath_chromatograms.getRunIDs() ) == 1
    rid = list(swath_chromatograms.getRunIDs())[0]

    start = time.time()
    initial_alignment_cutoff = 0.0001
    max_rt_diff = 30
    sd_data = -1 # We do not use the standard deviation data in this algorithm
    tr_data = transformations.LightTransformationData()
    spl_aligner = SplineAligner(initial_alignment_cutoff, new_exp)

    if method == "singleClosestRun":
        tree_mapped = None

        run_1 = [r for r in new_exp.runs if r.get_id() == rid][0]
        dist_matrix = getDistanceMatrix(new_exp, multipeptides, spl_aligner, singleRowId=run_1.get_id())
        print("Distance matrix took %ss" % (time.time() - start) )

        start = time.time()
        for run_0 in new_exp.runs:
            helper.addDataToTrafo(tr_data, run_0, run_1, spl_aligner, multipeptides,
                options.realign_method, max_rt_diff, sd_max_data_length=sd_data, force=is_test)

    elif method == "singleShortestPath":
        dist_matrix = None

        tree = MinimumSpanningTree(getDistanceMatrix(new_exp, multipeptides, spl_aligner))
        tree_mapped = [(new_exp.runs[a].get_id(), new_exp.runs[b].get_id()) for a,b in tree]
        print("Distance matrix took %ss" % (time.time() - start) )

        start = time.time()
        for edge in tree:
            helper.addDataToTrafo(tr_data, new_exp.runs[edge[0]], 
                new_exp.runs[edge[1]], spl_aligner, multipeptides, 
                options.realign_method, max_rt_diff, sd_max_data_length=sd_data, force=is_test)

    else:
        raise Exception("Unknown method: " + method)

    print("Alignment took %ss" % (time.time() - start) )
    start = time.time()
    multipeptides = analyze_multipeptides(new_exp, multipeptides, swath_chromatograms,
        tr_data, options.border_option, rid, tree=tree_mapped, mat=dist_matrix,
        disable_isotopic_transfer=options.disable_isotopic_transfer, is_test=is_test)
    print("Analyzing the runs took %ss" % (time.time() - start) )

    return new_exp, multipeptides, rid

def runImputeValues(options, peakgroups_file, trafo_fnames, is_test):
    """Impute values across chromatograms

    Args:
        peakgroups_file(filename): CSV file containing all peakgroups
        trafo_fnames(filename): A list of .tr filenames (it is assumed that in
            the same directory also the chromatogram mzML reside)
    Returns:
        A tuple of:
            new_exp(AlignmentExperiment): experiment containing the aligned peakgroups
            multipeptides(list(AlignmentHelper.Multipeptide)): list of multipeptides
            rid(string): run id of the analyzed run

    This function will read the csv file with all peakgroups as well as the
    transformation files (.tr) and the corresponding raw chromatograms which
    need to be in the same folder. It will then try to impute missing values
    for those peakgroups where no values is currently present, reading the raw
    chromatograms.
    """

    # We do not want to exclude any peakgroups for noiseIntegration (we assume
    # that alignment has already happened)
    fdr_cutoff_all_pg = 1.0

    start = time.time()
    reader = SWATHScoringReader.newReader([peakgroups_file],
                                          options.file_format,
                                          readmethod="complete",
                                          enable_isotopic_grouping = not options.disable_isotopic_grouping, 
                                          read_cluster_id=True)
    new_exp = Experiment()
    new_exp.runs = reader.parse_files()
    multipeptides = new_exp.get_all_multipeptides(fdr_cutoff_all_pg, verbose=False)
    print("Parsing the peakgroups file took %ss" % (time.time() - start) )

    start = time.time()
    transformation_collection_ = transformations.TransformationCollection()
    for filename in trafo_fnames:
        transformation_collection_.readTransformationData(filename)

    # Read the datapoints and perform the smoothing
    print("Reading the trafo file took %ss" % (time.time() - start) )
    start = time.time()
    transformation_collection_.initialize_from_data(reverse=True, smoother=options.realign_method, force=is_test)
    print("Initializing the trafo file took %ss" % (time.time() - start) )

    rid = None
    if options.do_single_run and not options.dry_run:
        # Do only a single run : read only one single file
        start = time.time()
        swath_chromatograms = SwathChromatogramCollection()
        swath_chromatograms.parseFromTrafoFiles([ options.do_single_run ])
        print("Reading the chromatogram files took %ss" % (time.time() - start) )
        assert len(swath_chromatograms.getRunIDs() ) == 1
        rid = swath_chromatograms.getRunIDs()[0]
        # 
        start = time.time()
        multipeptides = analyze_multipeptides(new_exp, multipeptides, swath_chromatograms,
            transformation_collection_, options.border_option, 
            onlyExtractFromRun=rid, disable_isotopic_transfer=options.disable_isotopic_transfer, is_test=is_test)
        print("Analyzing the runs took %ss" % (time.time() - start) )
        return new_exp, multipeptides, rid

    swath_chromatograms = SwathChromatogramCollection()
    swath_chromatograms.parseFromTrafoFiles(trafo_fnames)
    print("Reading the chromatogram files took %ss" % (time.time() - start) )

    if options.dry_run:
        print("Dry Run only")
        print("Found multipeptides:", len(multipeptides))
        print("Found swath chromatograms:", swath_chromatograms)
        return [], [], rid

    start = time.time()
    if options.cache_in_memory:
        run_ids = [r.get_id() for r in new_exp.runs]
        for rid in run_ids:
            # Create the cache for run "rid" and then only extract peakgroups from this run
            swath_chromatograms.createRunCache(rid)
            multipeptides = analyze_multipeptides(new_exp, multipeptides, 
                swath_chromatograms, transformation_collection_, options.border_option, 
                onlyExtractFromRun=rid, disable_isotopic_transfer=options.disable_isotopic_transfer, is_test=is_test)
    else:
        multipeptides = analyze_multipeptides(new_exp, multipeptides, swath_chromatograms,
            transformation_collection_, options.border_option, disable_isotopic_transfer=options.disable_isotopic_transfer, is_test=is_test)
    print("Analyzing the runs took %ss" % (time.time() - start) )
    return new_exp, multipeptides, rid

def analyze_multipeptides(new_exp, multipeptides, swath_chromatograms, 
                          transformation_collection_, border_option,
                          onlyExtractFromRun=None, tree=None, mat=None,
                          disable_isotopic_transfer=False, is_test=False):
    """Analyze the multipeptides and impute missing values

    This function has three different modes:

        - if a tree is given, it will use the single shortest path on the tree to infer the boundaries
        - if a distance matrix is given, it will use the closest run overall to infer the boundaries
        - if neither is given, it will use a mean / median / maximum of all
          borders by converting all boundaries to the frame of the reference run
          and then to the current run

    Args:
        new_exp(AlignmentExperiment): experiment containing the aligned peakgroups
        multipeptides(list(AlignmentHelper.Multipeptide)): list of multipeptides
        swath_chromatograms(dict): containing the objects pointing to the original chrom mzML (see runImputeValues)
        transformation_collection_(.TransformationCollection): specifying how to transform between retention times of different runs
        border_option(String): one of the following options ("mean", "median" "max_width"), determining how to aggregate multiple peak boundary information
        onlyExtractFromRun(String): Whether only to perform signal extraction from a single run (if not None, needs to be a run id)
        tree(MinimumSpanningTree): alignment guidance tree (if given, will use shortest path approach)
        mat( matrix(float) ): distance matrix, see getDistanceMatrix (if given, will use closest run approach)
        disable_isotopic_transfer(bool): whether to use isotopic grouping (e.g. group heavy/light channels together)

    Returns:
        The updated multipeptides

    This function will update the input multipeptides and add peakgroups, imputing missing values 
    """

    # Go through all aligned peptides
    class CounterClass(object): pass
    cnt = CounterClass()
    cnt.integration_bnd_warnings = 0
    cnt.imputations = 0
    cnt.imputation_succ = 0
    cnt.peakgroups = 0

    print("Will work on %s peptides" % (len(multipeptides)))
    for i,m in enumerate(multipeptides):
        if i % 500 == 0: 
            print("Done with %s out of %s" % (i, len(multipeptides)))
        if VERBOSE:
            print("Do group of peptides", m.getPrecursorGroups()[0].getPeptideGroupLabel())


        # Ensure that each run has either exactly one peakgroup (or zero):
        # If all peakgroups are old-style peakgroups (only single cluster) then
        # all peptides can only have a single peakgroup.
        if all([pg.get_cluster_id() == -1 for p in m.getAllPeptides() for pg in p.get_all_peakgroups()]):
            if any([len(p.peakgroups) > 1 for p in m.getAllPeptides()] ):
                raise Exception("Found more than one peakgroup for peptide %s - \
                    \n is this after alignment or did you forget to run feature_alignment.py beforehand?" % (
                    m.get_id() ))

        clusters = set( [pg.get_cluster_id() for p in m.getAllPeptides() for pg in p.get_all_peakgroups()])
        if VERBOSE: print("Shall work on clusters ", clusters)
        for cl in clusters:

            # Retrieve all transition group ids available for this precursor group
            # Iterate through all ids (in order of ascending FDR since we
            # believe transferring borders between runs is most efficient when
            # first using the group with the best overall score).
            trgr_ids = set([ trgr.get_id() for prgr in m.getPrecursorGroups() for trgr in prgr ])
            fdr_trgr_mapping = []
            for trgr_id in trgr_ids:
                selected_pg = [pg.get_fdr_score() for p in m.getAllPeptides() for pg in p.get_all_peakgroups()
                    if pg.peptide.get_id() == trgr_id and pg.get_cluster_id() == cl]

                # In some cases there is no peakgroup for the current cluster available
                if len(selected_pg) > 0:
                    fdr_trgr_mapping.append( (min(selected_pg), trgr_id) )

            if VERBOSE: print("mapping", sorted(fdr_trgr_mapping), "vs", trgr_ids)
            for fdr, trgr_id in sorted(fdr_trgr_mapping):
                if VERBOSE: print("================================================")
                if VERBOSE: print("Shall work on transition group", trgr_id, "with cluster", cl)

                # Get all selected peakgroups that correspond to the current transition group id 
                selected_pg = [pg for p in m.getAllPeptides() for pg in p.get_all_peakgroups()
                    if pg.peptide.get_id() == trgr_id and pg.get_cluster_id() == cl]
                analyze_multipeptide_cluster(m, cnt, new_exp, swath_chromatograms, 
                                          transformation_collection_, border_option, selected_pg, cl,
                                          onlyExtractFromRun, tree, mat, is_test=is_test)


            # Ensure consistency in peak picking between isotopic channels by
            # looking at the best channel for each run and then forcing the
            # peak boundaries of this channel unto the other channels. This is
            # done even if a very good peak was found in the other channels,
            # but ensuring that the peak boundaries are the same is more
            # important here.
            # disable_isotopic_transfer turns this feature off
            if disable_isotopic_transfer:
                continue

            for rid in [r.get_id() for r in new_exp.runs]:

                # Skip if we should not extract from this run
                if not onlyExtractFromRun is None:
                    if onlyExtractFromRun != rid:
                        continue

                if m.hasPrecursorGroup(rid):
                    prgr = m.getPrecursorGroup(rid)
                    pgs = [(pg_.get_fdr_score(), pg_) for prec_ in prgr for pg_ in prec_.peakgroups if pg_.get_cluster_id() == cl] 
                    if len(pgs) > 0:
                        best_pg = min(pgs)[1]
                        border_l = float(best_pg.get_value("leftWidth"))
                        border_r = float(best_pg.get_value("rightWidth"))
                        for prec in m.getPrecursorGroup(rid):
                            currpg = [pg for pg in prec.get_all_peakgroups() if pg.get_cluster_id() == cl]
                            if len(currpg) == 1 and currpg[0].get_feature_id() != best_pg.get_feature_id():
                                pg = currpg[0]
                                # Propagate peak boundaries of the best peakgroup to the current one
                                current_run = [r for r in new_exp.runs if r.get_id() == rid][0]
                                newpg = integrate_chromatogram(pg, current_run, swath_chromatograms,
                                                             border_l, border_r, cnt, is_test)
                                pg.set_value("Intensity", newpg.get_value("Intensity"))
                                pg.set_intensity(newpg.get_intensity())
                                pg.set_value("leftWidth", newpg.get_value("leftWidth"))
                                pg.set_value("rightWidth", newpg.get_value("rightWidth"))
                                pg.set_value("RT", newpg.get_value("RT"))
                                pg.set_value("m_score", 1.5)


        if all ([pg.get_cluster_id() == -1 for p in m.getAllPeptides() for pg in p.get_all_peakgroups()]):
            # no cluster info was read in -> select all
            for p in m.getAllPeptides():
                for pg in p.get_all_peakgroups():
                    pg.select_this_peakgroup()

    print("Imputations:", cnt.imputations, "Successful:", cnt.imputation_succ, "Still missing", cnt.imputations - cnt.imputation_succ)
    print("WARNING: %s times: Chromatogram does not cover full range"  % (cnt.integration_bnd_warnings))
    print("Peakgroups:", cnt.peakgroups)
    return multipeptides 

def analyze_multipeptide_cluster(current_mpep, cnt, new_exp, swath_chromatograms, 
                          transformation_collection_, border_option, selected_pg, cluster_id,
                          onlyExtractFromRun=None, tree=None, mat=None, is_test=False):

        for rid in [r.get_id() for r in new_exp.runs]:
            if VERBOSE: print("Shall work on run ", rid)
            cnt.peakgroups += 1

            # Check whether we have a peakgroup for this run id, if yes we mark
            # this peakgroup as selected for output, otherwise we try to impute.
            if not any([pg.peptide.run.get_id() == rid for pg in selected_pg]):

                # Skip if we should not extract from this run
                if not onlyExtractFromRun is None:
                    if onlyExtractFromRun != rid:
                        continue

                cnt.imputations += 1

                # Select current run, compute right/left integration border and then integrate
                current_run = [r for r in new_exp.runs if r.get_id() == rid][0]
                rmap = dict([(r.get_id(),i) for i,r in enumerate(new_exp.runs) ])

                if VERBOSE:
                    print("Will try to fill NA in run", current_run.get_id(), "for peptide", selected_pg[0].peptide.get_id())
                    print(tree, mat)

                inferred_from_isotope = False
                if current_mpep.hasPrecursorGroup(current_run.get_id()):
                    # If another precursor was already picked for this run,
                    # check if any peakgroups of the same cluster exist and if
                    # so, try to use the one with the minimal FDR score to
                    # infer the boundaries.
                    prgr = current_mpep.getPrecursorGroup(current_run.get_id())
                    pgs = [(pg_.get_fdr_score(), pg_) for prec_ in prgr for pg_ in prec_.peakgroups if pg_.get_cluster_id() == cluster_id] 
                    if len(pgs) > 0:
                        # TODO if the boundaries were inferred (e.g. fdr score > 1.0), should we not use them? 
                        best_pg = min(pgs)[1]
                        border_l = float(best_pg.get_value("leftWidth"))
                        border_r = float(best_pg.get_value("rightWidth"))
                        inferred_from_isotope = True

                        if VERBOSE:
                            print("Will try to infer from isotopically modified version of the same run!")
                            print("Precursor is", prgr)
                            for prec in prgr:
                                print(" --", prec)
                                for pg in prec.peakgroups:
                                    print("  * ", pg)
                            print("Min fdr pg", best_pg, " border %s / %s" % (border_l, border_r))

                if inferred_from_isotope:
                    # All good
                    pass

                elif tree is not None:
                    ## Use the closest path approach
                    border_l, border_r = integrationBorderShortestPath(selected_pg, 
                        rid, transformation_collection_, tree)
                elif mat is not None:
                    ## Use the closest overall run approach
                    border_l, border_r = integrationBorderShortestDistance(selected_pg, 
                        rid, transformation_collection_, mat, rmap)
                else:
                    ## Use the reference-based approach
                    border_l, border_r = integrationBorderReference(new_exp, selected_pg, rid, transformation_collection_, border_option)
                newpg = integrate_chromatogram(selected_pg[0], current_run, swath_chromatograms,
                                             border_l, border_r, cnt, is_test)
                if newpg != "NA": 
                    if VERBOSE: 
                        print("Managed to fill NA in run", current_run.get_id(), \
                          "with value", newpg.get_value("Intensity"), "/ borders", border_l, border_r)#, "for cluster", newpg.get_value("align_clusterid")
                    cnt.imputation_succ += 1
                    transition_group_id = newpg.get_value("transition_group_id")
                    try:
                        peptide_label = newpg.get_value("peptide_group_label")
                    except KeyError:
                        peptide_label = transition_group_id

                    # Create new precursor
                    precursor = GeneralPrecursor(transition_group_id, current_run)
                    newpg.setClusterID(cluster_id)
                    newpg.peptide = precursor
                    precursor.add_peakgroup(newpg)
                    precursor.sequence = selected_pg[0].peptide.sequence
                    precursor.protein_name = selected_pg[0].peptide.protein_name

                    if current_mpep.hasPrecursorGroup(current_run.get_id()):
                        prec_group = current_mpep.getPrecursorGroup(current_run.get_id())
                        if prec_group.getPrecursor(transition_group_id) is None:
                            # No precursors exists yet for this transition_group_id - this
                            # means that we have a new run for which a precursor group already
                            # exists (e.g. we have found a light version already) and are now
                            # dealing with the heavy ...
                            prec_group.addPrecursor(precursor)
                        else:
                            # Likely, another peakgroup from the same run and peptide was
                            # picked already but not from the same cluster
                            prec_group.getPrecursor(transition_group_id).add_peakgroup(newpg)
                    else:
                        # Create new precursor group and insert
                        precursor_group = PrecursorGroup(peptide_label, current_run)
                        precursor_group.addPrecursor(precursor)
                        current_mpep.insert(rid, precursor_group)

def integrate_chromatogram(template_pg, current_run, swath_chromatograms, 
                           left_start, right_end, cnt, is_test):
    """ Integrate a chromatogram from left_start to right_end and store the sum.

    Args:
        template_pg(GeneralPeakGroup): A template peakgroup from which to construct the new peakgroup
        current_run(SWATHScoringReader.Run): current run where the missing value occured
        swath_chromatograms(dict): containing the objects pointing to the original chrom mzML
        left_start(float): retention time for integration (left border)
        right_end(float): retention time for integration (right border)

    Returns:
        A new GeneralPeakGroup which contains the new, integrated intensity for this run (or "NA" if no chromatograms could be found).

    Create a new peakgroup from the old pg and then store the integrated intensity.
    """

    current_rid = current_run.get_id()
    orig_filename = current_run.get_openswath_filename()
    aligned_filename = current_run.get_aligned_filename()

    # Create new peakgroup by copying the old one
    newrow = ["NA" for ele in template_pg.row]
    newpg = GeneralPeakGroup(newrow, current_run, None)
    for element in ["transition_group_id"]:
        newpg.set_value(element, template_pg.get_value(element))
    for element in ["decoy", "Sequence", "FullPeptideName", "Charge", "ProteinName", "nr_peaks", "m.z", "m/z", "align_clusterid", "peptide_group_label"]:
        try:
            newpg.set_value(element, template_pg.get_value(element))
        except KeyError:
            # These are non-essential column names, we can just skip them ...
            pass

    newpg.set_value("align_runid", current_rid)
    newpg.set_value("run_id", current_rid)
    newpg.set_value("filename", orig_filename)
    newpg.set_value("align_origfilename", aligned_filename)
    newpg.set_value("RT", (left_start + right_end) / 2.0)
    newpg.set_value("leftWidth", left_start)
    newpg.set_value("rightWidth", right_end)
    newpg.set_value("m_score", 2.0)
    newpg.set_value("d_score", '-10') # -10 has a p value of 1.0 for 1-right side cdf
    newpg.set_value("peak_group_rank", -1)

    import uuid 
    thisid = str(uuid.uuid1())
    newpg.set_normalized_retentiontime((left_start + right_end) / 2.0)
    newpg.set_fdr_score(2.0)
    newpg.set_feature_id(thisid)
    if not is_test:
        newpg.set_value("id", thisid)

    integrated_sum = 0
    chrom_ids = template_pg.get_value("aggr_Fragment_Annotation").split(";")
    peak_areas = []

    # Sanity check: there should be the same number of fragment ion /
    #               chromatogram ids as there are peak areas. If this fails
    #               then the data is inconsistent and it is likely something
    #               went wrong, for example the fragment ids may have ";" in
    #               them.
    try:
        old_peak_areas = template_pg.get_value("aggr_Peak_Area").split(";")
        if len(old_peak_areas) != len(chrom_ids):
            print("WARNING: something is wrong, got %s fragment ids but %s fragment peak areas" % (len(chrom_ids), len(old_peak_areas)))
    except:
        # Caught the exception
        pass

    for chrom_id in chrom_ids:
        chromatogram = swath_chromatograms.getChromatogram(current_rid, chrom_id)
        # Check if we got a chromatogram with at least 1 peak: 
        if chromatogram is None:
            print("chromatogram is None (tried to get %s from run %s)" % (chrom_id, current_rid))
            return("NA")
        if len(chromatogram.peaks) == 0:
            print("chromatogram has no peaks None (tried to get %s from run %s)" % (chrom_id, current_rid))
            return("NA")

        # Check whether we even have data between left_start and right_end ... 
        if chromatogram.peaks[0][0] > left_start or chromatogram.peaks[-1][0] < right_end:
            cnt.integration_bnd_warnings += 1
            # print "WARNING: Chromatogram (%s,%s) does not cover full range (%s,%s)"  % (
            #     chromatogram.peaks[0][0],chromatogram.peaks[-1][0], left_start, right_end)

        curr_int = sum( [p[1] for p in chromatogram.peaks if p[0] > left_start and p[0] < right_end ])
        peak_areas.append(curr_int)
        # print integrated_sum, "chromatogram", sum(p[1] for p in chromatogram.peaks), \
        # "integrated from \t%s\t%s in run %s" %( left_start, right_end, current_rid), newpg.get_value("transition_group_id")

    integrated_sum = sum(peak_areas)
    aggr_annotation = ";".join(chrom_ids)
    aggr_peak_areas = ";".join([str(it) for it in peak_areas])
    newpg.set_value("aggr_Peak_Area", aggr_peak_areas)
    newpg.set_value("aggr_Fragment_Annotation", aggr_annotation)
    newpg.set_value("Intensity", integrated_sum)
    newpg.set_intensity(integrated_sum)
    return newpg

def write_out(new_exp, multipeptides, outfile, matrix_outfile, single_outfile=None):
    """ Write the result to disk 

    This writes all peakgroups to disk (newly imputed ones as previously found
    ones) as even some "previously good" peakgroups may have changed location
    due to isotopic_transfer.
    """
    # write out the complete original files 
    writer = csv.writer(open(outfile, "w"), delimiter="\t")
    header_first = new_exp.runs[0].header
    for run in new_exp.runs:
        assert header_first == run.header
    writer.writerow(header_first)

    print("number of precursors quantified:", len(multipeptides))
    for m in sorted(multipeptides, key=lambda x: str(x)):
        # selected_peakgroups = [p.peakgroups[0] for p in m.get_peptides()]
        # if (len(selected_peakgroups)*2.0 / len(new_exp.runs) < fraction_needed_selected) : continue
        for p in m.getAllPeptides():
            for selected_pg in sorted(p.peakgroups):
                if single_outfile is not None:
                    if single_outfile == selected_pg.get_value("run_id"):
                        # Only write the values for this run ... 
                        row_to_write = selected_pg.row
                        writer.writerow(row_to_write)
                else:
                    row_to_write = selected_pg.row
                    writer.writerow(row_to_write)

    if len(matrix_outfile) > 0:
        helper.write_out_matrix_file(matrix_outfile, new_exp.runs, multipeptides, 0.0, style=options.matrix_output_method)

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program will impute missing values"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", nargs = '+', required=False, help = 'A list of transformation files in the same folder as the .chrom.mzML files')
    parser.add_argument("--peakgroups_infile", dest="peakgroups_infile", required=True, help="Infile containing peakgroups (outfile from feature_alignment.py)")
    parser.add_argument("--out", dest="output", required=True, help="Output file with imputed values")
    parser.add_argument('--file_format', default='openswath', help="Which input file format is used (openswath, mprophet or peakview)")
    parser.add_argument("--out_matrix", dest="matrix_outfile", default="", help="Matrix containing one peak group per row (supports .csv, .tsv or .xls)")
    parser.add_argument("--matrix_output_method", dest="matrix_output_method", default='none', help="Which columns are written besides Intensity (none, RT, score, source or full)")
    parser.add_argument('--border_option', default='median', metavar="median", help="How to determine integration border (possible values: max_width, mean, median). Max width will use the maximal possible width (most conservative since it will overestimate the background signal).")
    parser.add_argument('--dry_run', action='store_true', default=False, help="Perform a dry run only")
    parser.add_argument('--test', dest="is_test", action='store_true', default=False, help="For running the tests (does not add a random id to the results)")
    parser.add_argument('--cache_in_memory', action='store_true', default=False, help="Cache data from a single run in memory")
    parser.add_argument('--method', dest='method', default="reference", help="Which method to use (singleShortestPath, singleClosestRun, reference)")
    parser.add_argument('--realign_runs', dest='realign_method', default="lowess", help="How to re-align runs in retention time ('diRT': use only deltaiRT from the input file, 'linear': perform a linear regression using best peakgroups, 'splineR': perform a spline fit using R, 'splineR_external': perform a spline fit using R (start an R process using the command line), 'splinePy' use Python native spline from scikits.datasmooth (slow!), 'lowess': use Robust locally weighted regression (lowess smoother), 'nonCVSpline, CVSpline': splines with and without cross-validation, 'Earth' : use Multivariate Adaptive Regression Splines using py-earth")
    parser.add_argument('--verbosity', default=0, type=int, help="Verbosity")
    parser.add_argument('--do_single_run', default='', metavar="", help="Only do a single run")

    experimental_parser = parser.add_argument_group('experimental options')
    experimental_parser.add_argument('--disable_isotopic_grouping', action='store_true', default=False, help="Disable grouping of isotopic variants by peptide_group_label, thus disabling matching of isotopic variants of the same peptide across channels. If turned off, each isotopic channel will be matched independently of the other. If enabled, the more certain identification will be used to infer the location of the peak in the other channel.")
    experimental_parser.add_argument('--disable_isotopic_transfer', action='store_true', default=False, help="Disable the transfer of isotopic boundaries in all cases. If enabled (default), the best (best score) isotopic channel dictates the peak boundaries and all other channels use those boundaries. This ensures consistency in peak picking in all cases.")

    args = parser.parse_args(sys.argv[1:])
    args.verbosity = int(args.verbosity)
    return args

def main(options):

    rid = None
    if options.method in ["singleShortestPath", "singleClosestRun"]:

        # Some parameter checking
        if options.do_single_run == "":
            raise Exception("Input mzML file (--do_single_run) cannot be empty when choosing a tree-based approacht")

        new_exp, multipeptides, rid = runSingleFileImputation(options,
                                                              options.peakgroups_infile,
                                                              options.do_single_run,
                                                              options.method,
                                                              options.is_test)
    else:
        if not options.method == "reference" and not options.method == "allTrafo":
            print ("Found unknown method '%s', assuming 'reference'." % options.method)

        new_exp, multipeptides, rid = runImputeValues(options,
                                                      options.peakgroups_infile,
                                                      options.infiles,
                                                      options.is_test)

    if options.dry_run:
        return

    write_out(new_exp, multipeptides, options.output, options.matrix_outfile, rid)

if __name__=="__main__":
    options = handle_args()
    main(options)
