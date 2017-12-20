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
import os, sys, csv
import numpy
import argparse
from msproteomicstoolslib.math.chauvenet import chauvenet
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import write_out_matrix_file

class Experiment(MRExperiment):
    """
    An Experiment is a container for multiple experimental runs - some of which may contain the same precursors.
    """

    def __init__(self):
        super(Experiment, self).__init__()
        self.transformation_collection = TransformationCollection()

    def get_max_pg(self):
      return len(self.runs)*len(self.union_transition_groups_set)

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
                          if prec.find_best_peptide_pg().peptide.get_decoy()])
        d.nr_targets = sum([len(prec.get_selected_peakgroups()) for prec in precursors_to_be_used 
                          if not prec.find_best_peptide_pg().peptide.get_decoy()])
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

    def print_stats(self, multipeptides, alignment, outlier_detection, fdr_cutoff, fraction_present, min_nrruns):
        nr_precursors_total = len(self.union_transition_groups_set)

        # Do statistics and print out
        in_all_runs_wo_align = len([1 for m in multipeptides if m.all_above_cutoff(fdr_cutoff) and not m.get_decoy()])
        proteins_in_all_runs_wo_align = len(set([m.find_best_peptide_pg().peptide.protein_name for m in multipeptides if m.all_above_cutoff(fdr_cutoff)]))
        proteins_in_all_runs_wo_align_target = len(set([m.find_best_peptide_pg().peptide.protein_name for m in multipeptides if m.all_above_cutoff(fdr_cutoff) and not m.find_best_peptide_pg().peptide.get_decoy()]))
        nr_all_proteins = len(set([m.find_best_peptide_pg().peptide.protein_name for m in multipeptides if not m.find_best_peptide_pg().peptide.get_decoy()]))
        nr_all_peptides = len(set([m.find_best_peptide_pg().peptide.sequence for m in multipeptides if not m.find_best_peptide_pg().peptide.get_decoy()]))
        peptides_in_all_runs_wo_align_target = len(set([m.find_best_peptide_pg().peptide.sequence for m in multipeptides if m.all_above_cutoff(fdr_cutoff) and not m.find_best_peptide_pg().peptide.get_decoy()]))

        print("Present targets in all runs", in_all_runs_wo_align)
        precursors_in_all_runs = [m for m in multipeptides if m.all_selected()]
        precursors_quantified = [m for m in multipeptides if len(m.get_selected_peakgroups()) > 0]

        # precursors_in_all_runs = [m for m in multipeptides if m.all_selected()]
        # precursors_in_all_runs = [m for m in multipeptides if m.all_selected()]
        nr_decoys = len([1 for prec in precursors_in_all_runs if prec.find_best_peptide_pg().peptide.get_decoy()])

        decoy_precursors = len([1 for m in multipeptides if len(m.get_selected_peakgroups()) > 0 and m.find_best_peptide_pg().peptide.get_decoy()])

        nr_peptides = len(set([prec.find_best_peptide_pg().peptide.sequence for prec in precursors_in_all_runs]))
        nr_proteins = len(set([prec.find_best_peptide_pg().peptide.protein_name for prec in precursors_in_all_runs]))
        nr_peptides_target = len(set([prec.find_best_peptide_pg().peptide.sequence for prec in precursors_in_all_runs if not prec.find_best_peptide_pg().peptide.get_decoy()]))
        nr_proteins_target = len(set([prec.find_best_peptide_pg().peptide.protein_name for prec in precursors_in_all_runs if not prec.find_best_peptide_pg().peptide.get_decoy()]))

        nr_precursors_to_quant = len(set([ prec for prec in precursors_quantified if not prec.find_best_peptide_pg().peptide.get_decoy()]))
        nr_proteins_to_quant = len(set([ prec.find_best_peptide_pg().peptide.protein_name for prec in precursors_quantified if not prec.find_best_peptide_pg().peptide.get_decoy()]))
        nr_peptides_to_quant = len(set([ prec.find_best_peptide_pg().peptide.sequence for prec in precursors_quantified if not prec.find_best_peptide_pg().peptide.get_decoy()]))

        nr_precursors_in_all = len([1 for m in multipeptides if m.all_selected() and not m.get_decoy()])
        max_pg = self.get_max_pg()
        dstats = self.estimate_real_fdr(multipeptides, fraction_present)
        dstats_all = self.estimate_real_fdr(multipeptides, 1.0)
        print("="*75)
        print("="*75)
        print("Total we have", len(self.runs), "runs with", len(self.union_transition_groups_set),\
                "peakgroups quantified in at least %s run(s) above FDR %0.4f %%" % (min_nrruns, fdr_cutoff*100) + ", " + \
                "giving maximally nr peakgroups", max_pg)
        print("We were able to quantify", alignment.nr_quantified, "/", max_pg, "peakgroups of which we aligned", \
                alignment.nr_aligned, "and changed order of", alignment.nr_changed, "and could not align", alignment.could_not_align)

        print("We were able to quantify %s / %s precursors in %s runs, and %s in all runs (up from %s before alignment)" % (
          nr_precursors_to_quant, nr_precursors_total, min_nrruns, nr_precursors_in_all, in_all_runs_wo_align))
        print("We were able to quantify %s / %s peptides in %s runs, and %s in all runs (up from %s before alignment)" % (
          nr_peptides_to_quant, nr_all_peptides, min_nrruns, nr_peptides_target, peptides_in_all_runs_wo_align_target))
        print("We were able to quantify %s / %s proteins in %s runs, and %s in all runs (up from %s before alignment)" % (
          nr_proteins_to_quant, nr_all_proteins, min_nrruns, nr_proteins_target, proteins_in_all_runs_wo_align_target))

        # print "quant proteins", nr_proteins_to_quant

        # Get decoy estimates
        if len(precursors_in_all_runs) > 0:
            print("Decoy percentage of peakgroups that are fully aligned %0.4f %% (%s out of %s) which roughly corresponds to a peakgroup FDR of %s %%" % (
                dstats_all.decoy_pcnt, dstats_all.nr_decoys, dstats_all.nr_decoys + dstats_all.nr_targets, dstats_all.est_real_fdr*100))

            print("Decoy percentage of peakgroups that are partially aligned %1.4f %% (%s out of %s) which roughly corresponds to a peakgroup FDR of %s %%" % (
                dstats.decoy_pcnt, dstats.nr_decoys, dstats.nr_decoys + dstats.nr_targets, dstats.est_real_fdr*100))

            print("There were", decoy_precursors, "decoy precursors identified out of", nr_precursors_to_quant + decoy_precursors, "precursors which is %0.4f %%" % (decoy_precursors *100.0 / (nr_precursors_to_quant + decoy_precursors)))


        if outlier_detection is not None: 
            print("Outliers:", outlier_detection.nr_outliers, "outliers in", len(multipeptides), "peptides or", outlier_detection.outlier_pg, "peakgroups out of", alignment.nr_quantified, "changed", outlier_detection.outliers_changed)

    def write_to_file(self, multipeptides, options):

        infiles = options.infiles
        outfile = options.outfile
        matrix_outfile = options.matrix_outfile
        matrix_excelfile = options.matix_excel
        yaml_outfile = options.yaml_outfile
        ids_outfile = options.ids_outfile
        fraction_needed_selected = options.min_frac_selected
        file_format = options.file_format

        selected_pgs = []
        for m in multipeptides:
            selected_peakgroups = m.get_selected_peakgroups()
            if (len(selected_peakgroups)*1.0 / len(self.runs) < fraction_needed_selected) : continue
            for p in m.getAllPeptides():
                selected_pg = p.get_selected_peakgroup()
                if selected_pg is None: continue
                selected_pgs.append(selected_pg)
        selected_ids_dict = dict( [ (pg.get_feature_id(), pg) for pg in selected_pgs] )

        if len(ids_outfile) > 0:
            fh = open(ids_outfile, "w")
            id_writer = csv.writer(fh, delimiter="\t")
            for pg in selected_pgs:
                id_writer.writerow([pg.get_feature_id()])
            fh.close()
            del id_writer

        if len(matrix_outfile) > 0:
            write_out_matrix_file(matrix_outfile, self.runs, multipeptides, fraction_needed_selected)


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
                if (len(selected_peakgroups)*1.0 / len(self.runs) < fraction_needed_selected) : continue
                for p in m.getAllPeptides():
                    selected_pg = p.get_selected_peakgroup()
                    if selected_pg is None: continue
                    row_to_write = selected_pg.row
                    row_to_write += [selected_pg.run.get_id(), selected_pg.run.orig_filename]
                    writer.writerow(row_to_write)
        elif len(outfile) > 0 and file_format == "openswath":
            # only in openswath we have the ID and can go back to the original file ... 
            # write out the complete original files 
            writer = csv.writer(open(outfile, "w"), delimiter="\t")
            header_first = self.runs[0].header
            for run in self.runs:
                assert header_first == run.header
            header_first += ["align_runid", "align_origfilename"]
            writer.writerow(header_first)

            for file_nr, f in enumerate(infiles):
              header_dict = {}
              reader = csv.reader(open(f), delimiter="\t")
              header = next(reader)
              for i,n in enumerate(header):
                header_dict[n] = i
              for row in reader:
                  f_id = row[ header_dict["id"]]
                  if f_id in selected_ids_dict:
                      # Check the "id" and "transition_group_id" field. 
                      # Unfortunately the id can be non-unique, there we check both.
                      trgroup_id = selected_ids_dict[f_id].peptide.get_id()
                      unique_peptide_id = row[ header_dict["transition_group_id"]]
                      if unique_peptide_id == trgroup_id:
                          row_to_write = row
                          row_to_write += [selected_ids_dict[f_id].peptide.run.get_id(), f]
                          writer.writerow(row_to_write)
 
        # Print out trafo data
        trafo_fnames = []
        for current_run in self.runs:
          current_id = current_run.get_id()
          ref_id = self.transformation_collection.getReferenceRunID() 
          filename = os.path.join(os.path.dirname(current_run.orig_filename), "transformation-%s-%s.tr" % (current_id, ref_id) )
          trafo_fnames.append(filename)
          self.transformation_collection.writeTransformationData(filename, current_id, ref_id)
          self.transformation_collection.readTransformationData(filename)

        if len(yaml_outfile) > 0:
            import yaml
            myYaml = {"RawData" : [], "PeakGroupData" : [ outfile ],
                      "ReferenceRun" : self.transformation_collection.getReferenceRunID() }
            for current_run in self.runs:
                current_id = current_run.get_id()
                ref_id = self.transformation_collection.getReferenceRunID() 
                filename = os.path.join(os.path.dirname(current_run.orig_filename), "transformation-%s-%s.tr" % (current_id, ref_id) )
                dirpath = os.path.realpath(os.path.dirname(current_run.orig_filename))
                this = {"id" : current_id, "directory" : dirpath, "trafo_file" : os.path.realpath(filename)}
                myYaml["RawData"].append(this)
            open(yaml_outfile, 'w').write(yaml.dump({"AlignedSwathRuns" : myYaml}))

        return trafo_fnames

# Detect outliers in "good" groups

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program will select all peakgroups below the FDR cutoff in all files and try to align them to each other."
    usage += "\nIf only one file is given, it will act as peakgroup selector (best by m_score)" + \
            "\nand will apply the provided FDR cutoff."

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", required=True, nargs = '+', help = 'A list of mProphet output files containing all peakgroups (use quotes around the filenames)')
    parser.add_argument("--out_matrix", dest="matrix_outfile", default="", help="Matrix containing one peak group per (supports .csv, .tsv or .xls)")

    parser.add_argument("--frac_selected", dest="min_frac_selected", default=0.0, type=float, help="Do not write peakgroup if selected in less than this fraction of runs (range 0 to 1)", metavar='0')
    parser.add_argument('--file_format', default='openswath', help="Which input file format is used (openswath or peakview)")
    parser.add_argument('--output_method', default='none', help="Which columns are written besides Intensity (none, RT, score or full)")
    parser.add_argument("--readmethod", dest="readmethod", default="minimal", help="Read full or minimal transition groups (minimal,full)")
    parser.add_argument('--remove_requant_values', action='store_true', default=False)
    parser.add_argument('--aligner_mscore_threshold', type=float, default=1.0, help="cutoff used at alignment, for coloring realigned values in blue")


    args = parser.parse_args(sys.argv[1:])

    if args.min_frac_selected < 0.0 or args.min_frac_selected > 1.0:
        raise Exception("Argument frac_selected needs to be a number between 0 and 1.0")

    return args

def fix_input_fnames(options, runs):
    """ Fix the input filenames 

    Replaces the run.orig_filename for each run with the filename that we think
    was most likely the original one where the run originated from. This value
    is taken from align_origfilename which probably contains the name of the
    original mProphet output file.
    """
    inputfile_mapping = {}
    aligned_run_id_name = "align_runid"
    filename_name = "align_origfilename"
    for file_nr, f in enumerate(options.infiles):
      stdout.flush()
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

      # Check if runs are already aligned (only one input file and correct header)
      already_aligned = (len(options.infiles) == 1 and aligned_run_id_name in header_dict)

      if not already_aligned:
          raise Exception("Can only complete data matrix generation on fully aligned runs")

      for this_row in reader:
          runid = this_row[header_dict[aligned_run_id_name]]
          filename = os.path.basename( this_row[header_dict[filename_name]] )
          if len(filename) > 0 and filename != "NA":
              inputfile_mapping[ runid ] = filename

    # Apply the fixed filenames
    for r in runs:
        r.orig_filename = inputfile_mapping[ r.get_id() ]

def main(options):
    import time

    # Read the files
    start = time.time()
    reader = SWATHScoringReader.newReader(options.infiles, options.file_format, options.readmethod)
    runs = reader.parse_files(True)
    # Create experiment
    this_exp = MRExperiment()
    this_exp.set_runs(runs)
    print("Reading the input files took %ss" % (time.time() - start) )

    # Fix input filenames
    fix_input_fnames(options, runs)

    # Map the precursors across multiple runs, determine the number of
    # precursors in all runs without alignment.
    start = time.time()
    multipeptides = this_exp.get_all_multipeptides(1.0, verbose=True)
    print("Mapping the precursors took %ss" % (time.time() - start) )

    for m in multipeptides:

        # Error handling if somehow more than one peakgroup was selected ... 
        for p in m.getAllPeptides():
            p._fixSelectedPGError(fixMethod="BestScore")

        if len(m.get_selected_peakgroups() ) > 0:
            continue 

        for p in m.getAllPeptides():
            if len(list(p.get_all_peakgroups())) != 1:
                print(p)
                print(dir(p))
                print(p.get_run_id())
                for pg in p.get_all_peakgroups():
                    print (pg.print_out())
                print (len(list(p.get_all_peakgroups())))

            assert len(list(p.get_all_peakgroups())) == 1
            for pg in p.get_all_peakgroups():
               pg.select_this_peakgroup()

    start = time.time()
    if len(options.matrix_outfile) > 0:
        write_out_matrix_file(options.matrix_outfile, this_exp.runs, multipeptides,
                              options.min_frac_selected, style=options.output_method, 
                              write_requant = not options.remove_requant_values, aligner_mscore_treshold=options.aligner_mscore_threshold)
    print("Writing output took %ss" % (time.time() - start) )

if __name__=="__main__":
    options = handle_args()
    main(options)
