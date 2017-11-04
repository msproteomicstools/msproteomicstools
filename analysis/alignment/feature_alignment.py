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

import os, sys, csv, time
import numpy
import argparse
from msproteomicstoolslib.math.chauvenet import chauvenet
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection, LightTransformationData
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment
from msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm import AlignmentAlgorithm
from msproteomicstoolslib.algorithms.alignment.AlignmentMST import getDistanceMatrix, TreeConsensusAlignment
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import write_out_matrix_file, addDataToTrafo
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
from msproteomicstoolslib.algorithms.alignment.FDRParameterEstimation import ParamEst
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree

class AlignmentStatistics(object):

    def __init__(self): 

        self.nr_runs = -1
        self.nr_aligned = 0
        self.nr_changed = 0
        self.nr_quantified = 0
        self.nr_removed = 0

        self.nr_good_peakgroups = 0
        self.nr_good_precursors = 0
        self.good_peptides = set([])
        self.good_proteins = set([])

        self.nr_quant_precursors = 0
        self.quant_peptides = set([])
        self.quant_proteins = set([])

    def count(self, multipeptides, fdr_cutoff, runs, skipDecoy=True):


        astats = self

        for m in multipeptides:

            if m.get_decoy() and skipDecoy:
                continue

            astats.nr_quantified += len(m.get_selected_peakgroups())

            # Count how many precursors / peptides / proteins fall below the threshold
            if m.find_best_peptide_pg().get_fdr_score() < fdr_cutoff:
                astats.nr_good_precursors += 1
                astats.good_peptides.update([m.getAllPeptides()[0].getSequence()])
                astats.good_proteins.update([m.getAllPeptides()[0].getProteinName()])

            # Count how many precursors / peptides / proteins were quantified
            if len(m.get_selected_peakgroups()) > 0:
                astats.nr_quant_precursors += 1
                astats.quant_peptides.update([m.getAllPeptides()[0].getSequence()])
                astats.quant_proteins.update([m.getAllPeptides()[0].getProteinName()])

            for p in m.getAllPeptides():

                # Count how many peakgroups simply fall below the threshold
                if p.get_best_peakgroup().get_fdr_score() < fdr_cutoff:
                    astats.nr_good_peakgroups += 1

                if p.get_selected_peakgroup() is not None:

                    # Number of peakgroups that are different from the original
                    if p.get_best_peakgroup().get_feature_id() != p.get_selected_peakgroup().get_feature_id() \
                       and p.get_selected_peakgroup().get_fdr_score() < fdr_cutoff:
                        astats.nr_changed += 1
                    # Number of peakgroups that were added
                    if p.get_best_peakgroup().get_fdr_score() > fdr_cutoff:
                        astats.nr_aligned += 1

                # Best peakgroup exists and is not selected
                elif p.get_best_peakgroup() is not None \
                  and p.get_best_peakgroup().get_fdr_score() < fdr_cutoff:
                    astats.nr_removed += 1

        self.max_pg = self.nr_good_precursors * len(runs)
        self.nr_runs = len(runs)

    def to_yaml(self):

        r = {}
        r["NumberRuns"] = self.nr_runs
        r["QuantifiedPrecursors"] = self.nr_good_precursors
        r["QuantifiedPeakgroups"] = self.nr_quantified
        r["MaximalPeakgroups"] = self.max_pg
        r["AlignedPeakgroups"] = self.nr_aligned
        r["ChangedPeakgroups"] = self.nr_changed
        r["NotAlignedPeakgroups"] = self.max_pg - self.nr_quantified
        r["RemovedPeakgroups"] = self.nr_removed
        r["AmbigousPeakgroups"] = self.nr_ambiguous
        r["MultipleSuitablePeakgroups"] = self.nr_multiple_align
        allr ={}
        allr["Precursors"] = {}
        allr["Peptides"] = {}
        allr["Proteins"] = {}

        allr["Precursors"]["Quant"] = self.nr_quant_precursors
        allr["Precursors"]["Total"] = self.nr_good_precursors
        allr["Precursors"]["AllRuns"] = self.nr_precursors_in_all
        allr["Precursors"]["AllRunsNoAlign"] = self.precursors_in_all_runs_wo_align
        allr["Peptides"]["Quant"] = len(self.quant_peptides)
        allr["Peptides"]["Total"] = len(self.good_peptides)
        allr["Peptides"]["AllRuns"] = self.nr_peptides_target
        allr["Peptides"]["AllRunsNoAlign"] = self.peptides_in_all_runs_wo_align_target

        allr["Proteins"]["Quant"] = len(self.quant_proteins)
        allr["Proteins"]["Total"] = len(self.good_proteins)
        allr["Proteins"]["AllRuns"] = self.nr_proteins_target
        allr["Proteins"]["AllRunsNoAlign"] = self.proteins_in_all_runs_wo_align_target

        r["AllRuns"] = allr

        return r

class Experiment(MRExperiment):
    """
    An Experiment is a container for multiple experimental runs - some of which may contain the same precursors.
    """

    def __init__(self):
        super(Experiment, self).__init__()
        self.transformation_collection = TransformationCollection()

    def estimate_real_fdr(self, multipeptides, fraction_needed_selected):
        class DecoyStats(object):
            def __init__(self):
                self.est_real_fdr = 0.0
                self.nr_decoys = 0
                self.nr_targets = 0
                self.decoy_pcnt = 0.0
                self.est_real_fdr = 0.0

        d = DecoyStats()
        precursors_to_be_used = [m for m in multipeptides if m.more_than_fraction_selected(fraction_needed_selected)]

        # count the decoys
        d.nr_decoys = sum([len(prec.get_selected_peakgroups()) for prec in precursors_to_be_used
                          if prec.find_best_peptide_pg().getPeptide().get_decoy()])
        d.nr_targets = sum([len(prec.get_selected_peakgroups()) for prec in precursors_to_be_used
                          if not prec.find_best_peptide_pg().getPeptide().get_decoy()])
        # estimate the real fdr by calculating the decoy ratio and dividing it
        # by the decoy ration obtained at @fdr_cutoff => which gives us the
        # decoy in/decrease realtive to fdr_cutoff. To calculate the absolute
        # value, we multiply by fdr_cutoff again (which was used to obtain the
        # original estimated decoy percentage).
        if self.estimated_decoy_pcnt is None: return d
        if (d.nr_targets + d.nr_decoys) == 0: return d
        d.decoy_pcnt = (d.nr_decoys * 100.0 / (d.nr_targets + d.nr_decoys) )
        d.est_real_fdr = d.decoy_pcnt / self.estimated_decoy_pcnt * self.initial_fdr_cutoff
        return d

    def print_stats(self, multipeptides, fdr_cutoff, fraction_present, min_nrruns):

        alignment = AlignmentStatistics()
        alignment.count(multipeptides, fdr_cutoff, self.runs)

        # Count presence in all runs (before alignment)
        mpep_target = [m for m in multipeptides if not m.get_decoy()]
        precursors_in_all_runs_wo_align_data = [m for m in mpep_target if m.all_above_cutoff(fdr_cutoff)]
        precursors_in_all_runs_wo_align = len(precursors_in_all_runs_wo_align_data)
        proteins_in_all_runs_wo_align_target = len(set([m.find_best_peptide_pg().getPeptide().getProteinName() for m in precursors_in_all_runs_wo_align_data]))
        peptides_in_all_runs_wo_align_target = len(set([m.find_best_peptide_pg().getPeptide().getSequence() for m in precursors_in_all_runs_wo_align_data]))

        # Count presence in all runs (before alignment)
        target_precursors_in_all_runs = [m for m in mpep_target if m.all_selected()]
        nr_peptides_target = len(set([prec.find_best_peptide_pg().getPeptide().getSequence() for prec in target_precursors_in_all_runs]))
        nr_proteins_target = len(set([prec.find_best_peptide_pg().getPeptide().getProteinName() for prec in target_precursors_in_all_runs]))

        nr_precursors_in_all = len(target_precursors_in_all_runs)
        max_pg = alignment.nr_good_precursors * len(self.runs)
        dstats = self.estimate_real_fdr(multipeptides, fraction_present)
        dstats_all = self.estimate_real_fdr(multipeptides, 1.0)

        # Get single/multiple hits stats
        from itertools import groupby
        precursors_quantified = [m for m in multipeptides if len(m.get_selected_peakgroups()) > 0]
        target_quant_protein_list = [ prec.find_best_peptide_pg().getPeptide().getProteinName() for prec in precursors_quantified 
                                     if not prec.find_best_peptide_pg().getPeptide().get_decoy()]
        target_quant_protein_list.sort()
        nr_sh_target_proteins = sum( [len(list(group)) == 1 for key, group in groupby(target_quant_protein_list)] )
        nr_mh_target_proteins = sum( [len(list(group)) > 1 for key, group in groupby(target_quant_protein_list)] )

        ### Store for later (yaml output)
        alignment.nr_ambiguous = self.nr_ambiguous
        alignment.nr_multiple_align = self.nr_multiple_align
        alignment.precursors_in_all_runs_wo_align = precursors_in_all_runs_wo_align
        alignment.peptides_in_all_runs_wo_align_target = peptides_in_all_runs_wo_align_target
        alignment.proteins_in_all_runs_wo_align_target = proteins_in_all_runs_wo_align_target

        alignment.nr_precursors_in_all = nr_precursors_in_all
        alignment.nr_peptides_target = nr_peptides_target
        alignment.nr_proteins_target = nr_proteins_target

        #
        ###########################################################################
        #
        print("="*75)
        print("="*75)
        print("Total we have", len(self.runs), "runs with", alignment.nr_good_precursors,
              "peakgroups quantified in at least %s run(s) below m_score (q-value) %0.4f %%" % (min_nrruns, fdr_cutoff*100) + ", " +
              "giving maximally nr peakgroups", max_pg)
        print("We were able to quantify", alignment.nr_quantified, "/", max_pg, "peakgroups of which we aligned",
              alignment.nr_aligned)
        print("  The order of", alignment.nr_changed, "peakgroups was changed,", max_pg - alignment.nr_quantified,
              "could not be aligned and %s were removed. Ambigous cases: %s, multiple suitable peakgroups: %s" % (
              alignment.nr_removed, self.nr_ambiguous, self.nr_multiple_align))
        print("We were able to quantify %s / %s precursors in %s runs, and %s in all runs (up from %s before alignment)" % (
          alignment.nr_quant_precursors, alignment.nr_good_precursors, min_nrruns, nr_precursors_in_all, precursors_in_all_runs_wo_align))
        print("We were able to quantify %s / %s peptides in %s runs, and %s in all runs (up from %s before alignment)" % (
              len(alignment.quant_peptides), len(alignment.good_peptides), min_nrruns, nr_peptides_target, peptides_in_all_runs_wo_align_target))
        print("We were able to quantify %s / %s proteins in %s runs, and %s in all runs (up from %s before alignment)" % (
              len(alignment.quant_proteins), len(alignment.good_proteins), min_nrruns, nr_proteins_target, proteins_in_all_runs_wo_align_target))
        print("Of these %s proteins, %s were multiple hits and %s were single hits." % (len(alignment.quant_proteins), nr_mh_target_proteins, nr_sh_target_proteins))

        # Get decoy estimates
        decoy_precursors = len([1 for m in multipeptides if len(m.get_selected_peakgroups()) > 0 and m.find_best_peptide_pg().getPeptide().get_decoy()])
        if len(target_precursors_in_all_runs) > 0:
            print("Decoy percentage of peakgroups that are fully aligned %0.4f %% (%s out of %s) which roughly corresponds to a peakgroup FDR of %s %%" % (
                dstats_all.decoy_pcnt, dstats_all.nr_decoys, dstats_all.nr_decoys + dstats_all.nr_targets, dstats_all.est_real_fdr*100))

            print("Decoy percentage of peakgroups that are partially aligned %0.4f %% (%s out of %s) which roughly corresponds to a peakgroup FDR of %s %%" % (
                dstats.decoy_pcnt, dstats.nr_decoys, dstats.nr_decoys + dstats.nr_targets, dstats.est_real_fdr*100))

            print("There were", decoy_precursors, "decoy precursors identified out of", \
                  alignment.nr_quant_precursors + decoy_precursors, "precursors which is %0.4f %%" % (
                      decoy_precursors * 100.0 / (alignment.nr_quant_precursors + decoy_precursors)))

        return alignment

    def _getTrafoFilename(self, current_run, ref_id):
        current_id = current_run.get_id()
        input_basename = os.path.basename(current_run.orig_filename)
        fn = os.path.splitext(input_basename)[0]
        dirname = os.path.dirname(current_run.orig_filename)
        filename = os.path.join(dirname, "%s-%s-%s.tr" % (fn, current_id, ref_id) )
        return filename

    def _write_trafo_files(self):
        # Print out trafo data
        trafo_fnames = []
        for current_run in self.runs:
            current_id = current_run.get_id()
            ref_id = self.transformation_collection.getReferenceRunID()
            filename = self._getTrafoFilename(current_run, ref_id)
            trafo_fnames.append(filename)
            self.transformation_collection.writeTransformationData(filename, current_id, ref_id)
            self.transformation_collection.readTransformationData(filename)

    def write_to_file(self, multipeptides, options, alignment, tree=None, writeTrafoFiles=True):

        infiles = options.infiles
        outfile = options.outfile
        matrix_outfile = options.matrix_outfile
        yaml_outfile = options.yaml_outfile
        ids_outfile = options.ids_outfile
        fraction_needed_selected = options.min_frac_selected
        file_format = options.file_format

        # 1. Collect ids of selected features
        selected_pgs = []
        for m in multipeptides:

            selected_peakgroups = m.get_selected_peakgroups()
            if (len(selected_peakgroups)*1.0 / len(self.runs)) < fraction_needed_selected: 
                continue

            for p in m.getAllPeptides():
                selected_pg = p.get_selected_peakgroup()
                clustered_pg = p.getClusteredPeakgroups()
                for pg in clustered_pg:
                    selected_pgs.append(pg)

        selected_ids_dict = dict( [ (pg.get_feature_id(), pg) for pg in selected_pgs] )

        # 2. Write out the (selected) ids
        if len(ids_outfile) > 0:
            fh = open(ids_outfile, "w")
            coll_ids = []
            for pg in selected_pgs:
                coll_ids.append(pg.get_feature_id())

            id_writer = csv.writer(fh, delimiter="\t")
            for pg in sorted(coll_ids):
                id_writer.writerow([pg])
            fh.close()
            del id_writer

        # 3. Write out the matrix outfile
        if len(matrix_outfile) > 0:
            write_out_matrix_file(matrix_outfile, self.runs, multipeptides,
                                  fraction_needed_selected,
                                  style=options.matrix_output_method,
                                  aligner_mscore_treshold=options.fdr_cutoff)

        # 4. Write out the full outfile
        if len(outfile) > 0 and options.readmethod == "full":
            # write out the complete original files 
            writer = csv.writer(open(outfile, "w"), delimiter="\t")
            header_first = self.runs[0].header
            for run in self.runs:
                assert header_first == run.header
            header_first += ["align_runid", "align_origfilename"]
            writer.writerow(header_first)

            for m in multipeptides:

                selected_peakgroups = m.get_selected_peakgroups()
                if (len(selected_peakgroups)*1.0 / len(self.runs)) < fraction_needed_selected:
                    continue

                for p in m.get_peptides():
                    selected_pg = p.get_selected_peakgroup()
                    if selected_pg is None: 
                        continue

                    row_to_write = selected_pg.row
                    row_to_write += [selected_pg.run.get_id(), selected_pg.run.orig_filename]
                    # Replace run_id with the aligned id (align_runid) ->
                    # otherwise the run_id is not guaranteed to be unique 
                    row_to_write[ header_dict["run_id"]] = selected_ids_dict[f_id].peptide.run.get_id()
                    writer.writerow(row_to_write)

        elif len(outfile) > 0 and file_format in ["openswath", "peakview_preprocess"]:

            name_of_id_col_map = { "openswath" : "id" , "peakview_preprocess" : "preprocess_id"}
            name_of_trgr_col_map = { "openswath" : "transition_group_id" , "peakview_preprocess" : "Pep Index"}
            name_of_id_col = name_of_id_col_map[file_format]
            name_of_trgr_col = name_of_trgr_col_map[file_format]

            # Only in openswath we have the ID and can go back to the original file.
            # We can write out the complete original files.

            writer = csv.writer(open(outfile, "w"), delimiter="\t")
            header_first = self.runs[0].header
            for run in self.runs:
                assert header_first == run.header
            header_first += ["align_runid", "align_origfilename", "align_clusterid"]
            writer.writerow(header_first)

            for file_nr, f in enumerate(infiles):
              header_dict = {}
              if f.endswith('.gz'):
                  import gzip
                  filehandler = gzip.open(f,'rb')
              else:
                  filehandler = open(f)

              reader = csv.reader(filehandler, delimiter="\t")
              header = next(reader)
              for i,n in enumerate(header):
                header_dict[n] = i

              for row in reader:
                  f_id = row[ header_dict[name_of_id_col]]
                  if f_id in selected_ids_dict:
                      # Check the "id" and "transition_group_id" field.
                      # Unfortunately the id can be non-unique, there we check both.
                      trgroup_id = selected_ids_dict[f_id].getPeptide().get_id()
                      unique_peptide_id = row[ header_dict[name_of_trgr_col]]
                      if unique_peptide_id == trgroup_id:
                          row_to_write = row
                          row_to_write += [selected_ids_dict[f_id].getPeptide().getRunId(), f, selected_ids_dict[f_id].get_cluster_id()]
                          # Replace run_id with the aligned id (align_runid) ->
                          # otherwise the run_id is not guaranteed to be unique 
                          if file_format == "openswath" : 
                              row_to_write[ header_dict["run_id"]] = selected_ids_dict[f_id].getPeptide().getRunId()
                          writer.writerow(row_to_write)

        # 5. Write out the .tr transformation files
        if writeTrafoFiles:
            self._write_trafo_files()

        # 6. Write out the YAML file
        if len(yaml_outfile) > 0:
            import yaml
            myYaml = {"Commandline" : sys.argv, 
                      "RawData" : [], "PeakGroupData" : [ outfile ],
                      "ReferenceRun" : self.transformation_collection.getReferenceRunID(), 
                      "FeatureAlignment" : 
                      {
                        "RawInputParameters" : options.__dict__,
                        "Parameters" : {}
                      },
                      "Parameters" : {}
                     }

            myYaml["Output"] = {}
            myYaml["Output"]["Tree"] = {}
            if tree is not None:
                myYaml["Output"]["Tree"]["Raw"] = [list(t) for t in tree]
                tree_mapped = [ [self.runs[a].get_id(), self.runs[b].get_id()] for a,b in tree]
                myYaml["Output"]["Tree"]["Mapped"] = tree_mapped
                tree_mapped = [ [self.runs[a].get_openswath_filename(), self.runs[b].get_openswath_filename()] for a,b in tree]
                myYaml["Output"]["Tree"]["MappedFile"] = tree_mapped
                tree_mapped = [ [self.runs[a].get_openswath_filename(), self.runs[b].get_openswath_filename()] for a,b in tree]
                myYaml["Output"]["Tree"]["MappedFile"] = tree_mapped
                tree_mapped = [ [self.runs[a].get_original_filename(), self.runs[b].get_original_filename()] for a,b in tree]
                myYaml["Output"]["Tree"]["MappedFileInput"] = tree_mapped

            myYaml["Output"]["Quantification"] = alignment.to_yaml()
            myYaml["Parameters"]["m_score_cutoff"] = float(options.fdr_cutoff) # deprecated
            myYaml["FeatureAlignment"]["Parameters"]["m_score_cutoff"] = float(options.fdr_cutoff)
            myYaml["FeatureAlignment"]["Parameters"]["fdr_cutoff"] = float(options.fdr_cutoff)
            myYaml["FeatureAlignment"]["Parameters"]["aligned_fdr_cutoff"] = float(options.aligned_fdr_cutoff)
            for current_run in self.runs:
                current_id = current_run.get_id()
                ref_id = self.transformation_collection.getReferenceRunID()
                filename = self._getTrafoFilename(current_run, ref_id)
                dirpath = os.path.dirname(current_run.orig_filename)
                ### Use real path (not very useful when moving data from one computer to another)
                ### filename = os.path.realpath(filename)
                ### dirpath = os.path.realpath(dirpath)
                this = {"id" : current_id, "directory" : dirpath, "trafo_file" : filename}
                myYaml["RawData"].append(this)
            open(yaml_outfile, 'w').write(yaml.dump({"AlignedSwathRuns" : myYaml}))

def estimate_aligned_fdr_cutoff(options, this_exp, multipeptides, fdr_range):
    print("Try to find parameters for target fdr %0.2f %%" % (options.target_fdr * 100))

    for aligned_fdr_cutoff in fdr_range:
        # Do the alignment and annotate chromatograms without identified features 
        # Then perform an outlier detection over multiple runs
        # unselect all
        for m in multipeptides:
            for p in m.get_peptides():
                p.unselect_all()

        # now align
        options.aligned_fdr_cutoff = aligned_fdr_cutoff
        alignment = align_features(multipeptides, options.rt_diff_cutoff, options.fdr_cutoff, options.aligned_fdr_cutoff, options.method)
        est_fdr = this_exp.estimate_real_fdr(multipeptides, options.min_frac_selected).est_real_fdr

        print("Estimated FDR: %0.4f %%" % (est_fdr * 100), "at position aligned fdr cutoff ", aligned_fdr_cutoff)
        if est_fdr > options.target_fdr:
            # Unselect the peptides again ...
            for m in multipeptides:
                for p in m.get_peptides():
                    p.unselect_all()
            return aligned_fdr_cutoff

def doMSTAlignment(exp, multipeptides, max_rt_diff, rt_diff_isotope, initial_alignment_cutoff,
                   fdr_cutoff, aligned_fdr_cutoff, smoothing_method, method,
                   use_RT_correction, stdev_max_rt_per_run, use_local_stdev, mst_use_ref, force, optimized_cython):
    """
    Minimum Spanning Tree (MST) based local aligment 
    """

    spl_aligner = SplineAligner(initial_alignment_cutoff, experiment=exp)

    if mst_use_ref:
        # force reference-based alignment
        bestrun = spl_aligner._determine_best_run(exp)
        ref = spl_aligner._determine_best_run(exp).get_id()
        refrun_id, refrun = [ (i,run) for i, run in enumerate(exp.runs) if run.get_id() == ref][0]
        tree = [( i, refrun_id) for i in range(len(exp.runs)) if i != refrun_id]
    else:
        start = time.time()
        tree = MinimumSpanningTree(getDistanceMatrix(exp, multipeptides, spl_aligner))
        print("Computing tree took %0.2fs" % (time.time() - start) )

    print("Computed Tree:", tree)

    
    # Get alignments
    start = time.time()
    try:
        from msproteomicstoolslib.cython._optimized import CyLightTransformationData
        if optimized_cython:
            tr_data = CyLightTransformationData()
        else:
            tr_data = LightTransformationData()
    except ImportError:
        print("WARNING: cannot import CyLightTransformationData, will use Python version (slower).")
        tr_data = LightTransformationData()

    for edge in tree:
        addDataToTrafo(tr_data, exp.runs[edge[0]], exp.runs[edge[1]],
                       spl_aligner, multipeptides, smoothing_method,
                       max_rt_diff, force=force)

    tree_mapped = [ (exp.runs[a].get_id(), exp.runs[b].get_id()) for a,b in tree]
    print("Computing transformations for all edges took %0.2fs" % (time.time() - start) )

    # Perform work
    al = TreeConsensusAlignment(max_rt_diff, fdr_cutoff, aligned_fdr_cutoff, 
                                rt_diff_isotope=rt_diff_isotope,
                                correctRT_using_pg=use_RT_correction,
                                stdev_max_rt_per_run=stdev_max_rt_per_run,
                                use_local_stdev=use_local_stdev)

    if method == "LocalMST":
        if optimized_cython:
            al.alignBestCluster(multipeptides, tree_mapped, tr_data)
        else:
            print("WARNING: cannot utilize optimized MST alignment (needs readmethod = cminimal), will use Python version (slower).")
            al.alignBestCluster_legacy(multipeptides, tree_mapped, tr_data)
    elif method == "LocalMSTAllCluster":
        al.alignAllCluster(multipeptides, tree_mapped, tr_data)

    # Store number of ambigous cases (e.g. where more than one peakgroup below
    # the strict quality cutoff was found in the RT window) and the number of
    # cases where multiple possibilities were found.
    exp.nr_ambiguous = al.nr_ambiguous
    exp.nr_multiple_align = al.nr_multiple_align

    return tree

def doParameterEstimation(options, this_exp, multipeptides):
    """
    Perform (q-value) parameter estimation
    """

    start = time.time()
    print("-"*35)
    print("Do Parameter estimation")
    p = ParamEst(min_runs=options.nr_high_conf_exp,verbose=True)
    decoy_frac = p.compute_decoy_frac(multipeptides, options.target_fdr)
    print("Found target decoy fraction overall %0.4f%%" % (decoy_frac*100))
    try:
        fdr_cutoff_calculated = p.find_iterate_fdr(multipeptides, decoy_frac)
    except UnboundLocalError:
        raise Exception("Could not estimate FDR accurately!")

    if fdr_cutoff_calculated > options.target_fdr:
        # Re-read the multipeptides with the new cutoff (since the cutoff might
        # be higher than before, we might have to consider more peptides now).
        multipeptides = this_exp.get_all_multipeptides(fdr_cutoff_calculated, verbose=True)
        print("Re-parse the files!")
        try:
            fdr_cutoff_calculated = p.find_iterate_fdr(multipeptides, decoy_frac)
        except UnboundLocalError:
            raise Exception("Could not estimate FDR accurately!")

    options.aligned_fdr_cutoff = float(options.aligned_fdr_cutoff)
    if options.aligned_fdr_cutoff < 0:
        # Estimate the aligned_fdr parameter -> if the new fdr cutoff is
        # lower than the target fdr, we can use the target fdr as aligned
        # cutoff but if its higher we have to guess (here we take
        # 2xcutoff).
        if fdr_cutoff_calculated < options.target_fdr:
            options.aligned_fdr_cutoff = options.target_fdr
        else:
            options.aligned_fdr_cutoff = 2*fdr_cutoff_calculated

    options.fdr_cutoff = fdr_cutoff_calculated
    print("Using an m_score (q-value) cutoff of %0.7f%%" % (fdr_cutoff_calculated*100))
    print("For the aligned values, use a cutoff of %0.7f%%" % (options.aligned_fdr_cutoff*100))
    print("Parameter estimation took %0.2fs" % (time.time() - start) )
    print("-"*35)
    return multipeptides

def doReferenceAlignment(options, this_exp, multipeptides):

    # Performing re-alignment using a reference run
    if options.realign_method != "diRT":
        start = time.time()
        spl_aligner = SplineAligner(alignment_fdr_threshold = options.alignment_score, 
                                   smoother=options.realign_method,
                                   external_r_tmpdir = options.tmpdir, 
                                   experiment=this_exp)
        this_exp.transformation_collection = spl_aligner.rt_align_all_runs(this_exp, multipeptides)
        trafoError = spl_aligner.getTransformationError()
        print("Aligning the runs took %0.2fs" % (time.time() - start) )

    try:
        options.aligned_fdr_cutoff = float(options.aligned_fdr_cutoff)
    except ValueError:
        # We have a range of values to step through. 
        # Since we trust the input, wo dont do error checking.
        exec("fdr_range = numpy.arange(%s)" % options.aligned_fdr_cutoff)
        options.aligned_fdr_cutoff = estimate_aligned_fdr_cutoff(options, this_exp, multipeptides, fdr_range)

    try:
        options.rt_diff_cutoff = float(options.rt_diff_cutoff)
    except ValueError:
        if options.rt_diff_cutoff == "auto_2medianstdev":
            options.rt_diff_cutoff = 2*numpy.median(list(trafoError.getStdev()))
        elif options.rt_diff_cutoff == "auto_3medianstdev":
            options.rt_diff_cutoff = 3*numpy.median(list(trafoError.getStdev()))
        elif options.rt_diff_cutoff == "auto_4medianstdev":
            options.rt_diff_cutoff = 4*numpy.median(list(trafoError.getStdev()))
        elif options.rt_diff_cutoff == "auto_maxstdev":
            options.rt_diff_cutoff = max(list(trafoError.getStdev()))
        else:
            raise Exception("max_rt_diff either needs to be a value in seconds or" + \
                            "one of ('auto_2medianstdev', 'auto_3medianstdev', " + \
                            "'auto_4medianstdev', 'auto_maxstdev'). Found instead: '%s'" % options.rt_diff_cutoff)

    print("Will calculate with aligned_fdr cutoff of", options.aligned_fdr_cutoff, "and an RT difference of", options.rt_diff_cutoff)
    start = time.time()
    AlignmentAlgorithm().align_features(multipeptides, 
                    options.rt_diff_cutoff, options.fdr_cutoff,
                    options.aligned_fdr_cutoff, options.method)
    print("Re-aligning peak groups took %0.2fs" % (time.time() - start) )

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program will select all peakgroups below the FDR cutoff in all files and try to align them to each other."
    usage += "\nIf only one file is given, it will act as peakgroup selector (best by m_score)" + \
            "\nand will apply the provided FDR cutoff."

    import ast
    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", required=True, nargs = '+', help = 'A list of mProphet output files containing all peakgroups (use quotes around the filenames)', metavar="INP")
    parser.add_argument('--file_format', default='openswath', help="Input file format (openswath, mprophet or peakview)", metavar="")
    parser.add_argument("--out", dest="outfile", required=True, default="feature_alignment_outfile", help="Output file with filtered peakgroups")
    parser.add_argument("--out_matrix", dest="matrix_outfile", default="", help="Matrix containing one peak group per row (supports .csv, .tsv or .xlsx)", metavar="")
    parser.add_argument("--out_ids", dest="ids_outfile", default="", help="Id file only containing the ids", metavar="")
    parser.add_argument("--out_meta", dest="yaml_outfile", default="", help="Outfile containing meta information", metavar="")
    parser.add_argument("--fdr_cutoff", dest="fdr_cutoff", default=0.01, type=float, help="Fixed FDR cutoff used for seeding (only assays where at least one peakgroup in one run is below this cutoff will be included in the result), see also target_fdr for a non-fixed cutoff", metavar='0.01')
    parser.add_argument("--target_fdr", dest="target_fdr", default=-1, type=float, help="If parameter estimation is used, which target FDR should be optimized for. If set to lower than 0, parameter estimation is turned off.", metavar='0.01')
    parser.add_argument("--max_fdr_quality", dest="aligned_fdr_cutoff", default=-1.0, help="Extension m-score score cutoff, peakgroups of this quality will still be considered for alignment during extension", metavar='-1')
    parser.add_argument("--max_rt_diff", dest="rt_diff_cutoff", default=30, help="Maximal difference in RT (in seconds) for two aligned features", metavar='30')
    parser.add_argument("--iso_max_rt_diff", dest="rt_diff_isotope", default=10, help="Maximal difference in RT (in seconds) for two isotopic channels in the same run", metavar='30')
    parser.add_argument("--frac_selected", dest="min_frac_selected", default=0.0, type=float, help="Do not write peakgroup if selected in less than this fraction of runs (range 0 to 1)", metavar='0')
    parser.add_argument('--method', default='best_overall', help="Method to use for the clustering (best_overall, best_cluster_score or global_best_cluster_score, global_best_overall, LocalMST, LocalMSTAllCluster).")
    parser.add_argument("--verbosity", default=0, type=int, help="Verbosity (0 = little)", metavar='0')
    parser.add_argument("--matrix_output_method", dest="matrix_output_method", default='none', help="Which columns are written besides Intensity (none, RT, score, source or full)", metavar="")
    parser.add_argument('--realign_method', dest='realign_method', default="diRT", help="RT alignment method (diRT, linear, splineR, splineR_external, splinePy, lowess, lowess_biostats, lowess_statsmodels, lowess_cython, nonCVSpline, CVSpline, Earth, WeightedNearestNeighbour, SmoothLLDMedian, None)", metavar="diRT")
    parser.add_argument('--force', action='store_true', default=False, help="Force alignment")

    mst_parser = parser.add_argument_group('options for the MST')

    mst_parser.add_argument("--mst:useRTCorrection", dest="mst_correct_rt", type=ast.literal_eval, default=True, help="Use aligned peakgroup RT to continue threading", metavar='True')
    mst_parser.add_argument("--mst:Stdev_multiplier", dest="mst_stdev_max_per_run", type=float, default=-1.0, help="How many standard deviations the peakgroup can deviate in RT during the alignment (if less than max_rt_diff, then max_rt_diff is used)", metavar='-1.0')
    mst_parser.add_argument("--mst:useLocalStdev", dest="mst_local_stdev", type=ast.literal_eval, default=False, help="Use local standard deviation of the  alignment", metavar='False')
    mst_parser.add_argument("--mst:useReference", dest="mst_use_ref", type=ast.literal_eval, default=False, help="Use a reference-based tree for alignment", metavar='False')

    experimental_parser = parser.add_argument_group('experimental options')

    experimental_parser.add_argument('--disable_isotopic_grouping', action='store_true', default=False, help="Disable grouping of isotopic variants by peptide_group_label")
    experimental_parser.add_argument('--use_dscore_filter', action='store_true', default=False)
    experimental_parser.add_argument("--dscore_cutoff", default=1.96, type=float, help="Discard all peakgroups below this d-score", metavar='1.96')
    experimental_parser.add_argument("--nr_high_conf_exp", default=1, type=int, help="Number of experiments in which the peptide needs to be identified with confidence above fdr_cutoff", metavar='1')
    experimental_parser.add_argument("--readmethod", dest="readmethod", default="minimal", help="Read full or minimal transition groups (cminimal,minimal,full)", metavar="minimal")
    experimental_parser.add_argument("--tmpdir", dest="tmpdir", default="/tmp/", help="Temporary directory")
    experimental_parser.add_argument("--alignment_score", dest="alignment_score", default=0.0001, type=float, help="Minimal score needed for a feature to be considered for alignment between runs", metavar='0.0001')

    # deprecated methods
    experimental_parser.add_argument('--realign_runs', action='store_true', default=False, help="Deprecated option (equals '--realign_method external_r')")
    experimental_parser.add_argument('--use_external_r', action='store_true', default=False, help="Deprecated option (equals '--realign_method external_r')")

    args = parser.parse_args(sys.argv[1:])

    # deprecated
    if args.realign_runs or args.use_external_r:
        print("WARNING, deprecated --realign_runs or --use_external_r used! Please use --realign_method instead")
        args.realign_method = "splineR_external"

    if args.min_frac_selected < 0.0 or args.min_frac_selected > 1.0:
        raise Exception("Argument frac_selected needs to be a number between 0 and 1.0")

    if args.target_fdr > 0:
        # Parameter estimation turned on: check user input ...
        if args.fdr_cutoff != 0.01:
            raise Exception("You selected parameter estimation with target_fdr - cannot set fdr_cutoff as well! It does not make sense to ask for estimation of the fdr_cutoff (target_fdr > 0.0) and at the same time specify a certain fdr_cutoff.")
        args.fdr_cutoff = args.target_fdr
        # if args.aligned_fdr_cutoff != -1.0:
        #     raise Exception("You selected parameter estimation with target_fdr - cannot set max_fdr_quality as well!")
        pass
    else:
        # Parameter estimation turned off: Check max fdr quality ...
        try:
            if float(args.aligned_fdr_cutoff) < 0:
                args.aligned_fdr_cutoff = args.fdr_cutoff
                print("Setting max_fdr_quality automatically to fdr_cutoff of", args.fdr_cutoff)
            elif float(args.aligned_fdr_cutoff) < args.fdr_cutoff:
                raise Exception("max_fdr_quality cannot be smaller than fdr_cutoff!")
        except ValueError:
            pass
    return args

def main(options):

    class DReadFilter(object):
        def __init__(self, cutoff):
            self.cutoff = cutoff
        def __call__(self, row, header):
            return float(row[ header["d_score" ] ]) > self.cutoff


    readfilter = ReadFilter()
    if options.use_dscore_filter:
        readfilter = DReadFilter(float(options.dscore_cutoff))

    # Read the files
    start = time.time()

    optimized_cython = options.realign_method in [ "splineR_external", "lowess", "lowess_biostats", "lowess_statsmodels", "lowess_cython"]
    optimized_cython = optimized_cython and options.readmethod == "cminimal"
    if optimized_cython:
        print("Provided arguments that allow execution of optimized cython code")

    reader = SWATHScoringReader.newReader(options.infiles, options.file_format,
                                          options.readmethod, readfilter,
                                          enable_isotopic_grouping = not options.disable_isotopic_grouping, 
                                          read_cluster_id = False)
    runs = reader.parse_files(options.realign_method != "diRT", options.verbosity, optimized_cython)

    # Create experiment
    this_exp = Experiment()
    this_exp.set_runs(runs)
    print("Reading the input files took %0.2fs" % (time.time() - start) )

    # Map the precursors across multiple runs, determine the number of
    # precursors in all runs without alignment.
    start = time.time()
    multipeptides = this_exp.get_all_multipeptides(options.fdr_cutoff, verbose=False, verbosity=options.verbosity)
    print("Mapping the precursors took %0.2fs" % (time.time() - start) )

    if options.target_fdr > 0:
        multipeptides = doParameterEstimation(options, this_exp, multipeptides)

    tree_out = None
    if options.method == "LocalMST" or options.method == "LocalMSTAllCluster":
        start = time.time()
        if options.mst_stdev_max_per_run > 0:
            stdev_max_rt_per_run = options.mst_stdev_max_per_run
        else:
            stdev_max_rt_per_run = None

        tree_out = doMSTAlignment(this_exp,
                       multipeptides, float(options.rt_diff_cutoff), 
                       float(options.rt_diff_isotope),
                       float(options.alignment_score), options.fdr_cutoff,
                       float(options.aligned_fdr_cutoff),
                       options.realign_method, options.method,
                       options.mst_correct_rt, stdev_max_rt_per_run,
                       options.mst_local_stdev, options.mst_use_ref, options.force, 
                       optimized_cython)
        print("Re-aligning peak groups took %0.2fs" % (time.time() - start) )
    else:
        doReferenceAlignment(options, this_exp, multipeptides)

    # Filter by high confidence (e.g. keep only those where enough high confidence IDs are present)
    start = time.time()
    for mpep in multipeptides:
        # check if we have found enough peakgroups which are below the cutoff
        count = 0
        for pg in mpep.get_selected_peakgroups():
            if pg.get_fdr_score() < options.fdr_cutoff:
                count += 1
        if count < options.nr_high_conf_exp:
            for p in mpep.getAllPeptides():
                p.unselect_all()
    print("Filtering took %0.2fs" % (time.time() - start) )

    # print statistics, write output
    start = time.time()
    al = this_exp.print_stats(multipeptides, options.fdr_cutoff, options.min_frac_selected, options.nr_high_conf_exp)
    this_exp.write_to_file(multipeptides, options, alignment=al, tree=tree_out)
    print("Writing output took %0.2fs" % (time.time() - start) )

if __name__=="__main__":
    options = handle_args()
    main(options)
