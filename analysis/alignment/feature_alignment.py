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

import os, sys, csv
import numpy
import argparse
from msproteomicstoolslib.math.chauvenet import chauvenet
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment
from msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm import AlignmentAlgorithm
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import write_out_matrix_file
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner

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
        class DecoyStats():
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

        print "Present targets in all runs", in_all_runs_wo_align
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
        nr_peptides_to_quant = len(set([ prec.find_best_peptide_pg().peptide.sequence for prec in precursors_quantified if not prec.find_best_peptide_pg().peptide.get_decoy()]))
        target_quant_protein_list = [ prec.find_best_peptide_pg().peptide.protein_name for prec in precursors_quantified if not prec.find_best_peptide_pg().peptide.get_decoy()]
        nr_proteins_to_quant = len(set(target_quant_protein_list))

        # TODO
        from itertools import groupby
        target_quant_protein_list.sort()
        nr_sh_target_proteins = sum( [len(list(group)) == 1 for key, group in groupby(target_quant_protein_list)] )
        nr_mh_target_proteins = sum( [len(list(group)) > 1 for key, group in groupby(target_quant_protein_list)] )

        nr_precursors_in_all = len([1 for m in multipeptides if m.all_selected() and not m.get_decoy()])
        max_pg = self.get_max_pg()
        dstats = self.estimate_real_fdr(multipeptides, fraction_present)
        dstats_all = self.estimate_real_fdr(multipeptides, 1.0)
        print "="*75
        print "="*75
        print "Total we have", len(self.runs), "runs with", len(self.union_transition_groups_set),\
                "peakgroups quantified in at least %s run(s) above FDR %0.4f %%" % (min_nrruns, fdr_cutoff*100) + ", " + \
                "giving maximally nr peakgroups", max_pg
        print "We were able to quantify", alignment.nr_quantified, "/", max_pg, "peakgroups of which we aligned", \
                alignment.nr_aligned, "and changed order of", alignment.nr_changed, "and could not align", alignment.could_not_align

        print "We were able to quantify %s / %s precursors in %s runs, and %s in all runs (up from %s before alignment)" % (
          nr_precursors_to_quant, nr_precursors_total, min_nrruns, nr_precursors_in_all, in_all_runs_wo_align)
        print "We were able to quantify %s / %s peptides in %s runs, and %s in all runs (up from %s before alignment)" % (
          nr_peptides_to_quant, nr_all_peptides, min_nrruns, nr_peptides_target, peptides_in_all_runs_wo_align_target)
        print "We were able to quantify %s / %s proteins in %s runs, and %s in all runs (up from %s before alignment)" % (
          nr_proteins_to_quant, nr_all_proteins, min_nrruns, nr_proteins_target, proteins_in_all_runs_wo_align_target)
        print "  Of these %s proteins, %s were multiple hits and %s were single hits" % (nr_proteins_to_quant, nr_mh_target_proteins, nr_sh_target_proteins)

        # print "quant proteins", nr_proteins_to_quant

        # Get decoy estimates
        if len(precursors_in_all_runs) > 0:
            print "Decoy percentage of peakgroups that are fully aligned %0.4f %% (%s out of %s) which roughly corresponds to a peakgroup FDR of %s %%" % (
                dstats_all.decoy_pcnt, dstats_all.nr_decoys, dstats_all.nr_decoys + dstats_all.nr_targets, dstats_all.est_real_fdr*100)

            print "Decoy percentage of peakgroups that are partially aligned %0.4f %% (%s out of %s) which roughly corresponds to a peakgroup FDR of %s %%" % (
                dstats.decoy_pcnt, dstats.nr_decoys, dstats.nr_decoys + dstats.nr_targets, dstats.est_real_fdr*100)

            print "There were", decoy_precursors, "decoy precursors identified out of", nr_precursors_to_quant + decoy_precursors, "precursors which is %0.4f %%" % (decoy_precursors *100.0 / (nr_precursors_to_quant + decoy_precursors))


        if outlier_detection is not None:
            print "Outliers:", outlier_detection.nr_outliers, "outliers in", len(multipeptides), "peptides or", outlier_detection.outlier_pg, "peakgroups out of", alignment.nr_quantified, "changed", outlier_detection.outliers_changed

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

    def write_to_file(self, multipeptides, options):

        infiles = options.infiles
        outfile = options.outfile
        matrix_outfile = options.matrix_outfile
        yaml_outfile = options.yaml_outfile
        ids_outfile = options.ids_outfile
        fraction_needed_selected = options.min_frac_selected
        file_format = options.file_format

        selected_pgs = []
        for m in multipeptides:
            selected_peakgroups = m.get_selected_peakgroups()
            if (len(selected_peakgroups)*1.0 / len(self.runs) < fraction_needed_selected) : continue
            for p in m.get_peptides():
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
                for p in m.get_peptides():
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
              if f.endswith('.gz'):
                  import gzip
                  filehandler = gzip.open(f,'rb')
              else:
                  filehandler = open(f)
              reader = csv.reader(filehandler, delimiter="\t")
              header = reader.next()
              for i,n in enumerate(header):
                header_dict[n] = i
              for row in reader:
                  f_id = row[ header_dict["id"]]
                  if selected_ids_dict.has_key(f_id):
                      # Check the "id" and "transition_group_id" field.
                      # Unfortunately the id can be non-unique, there we check both.
                      trgroup_id = selected_ids_dict[f_id].peptide.get_id()
                      unique_peptide_id = row[ header_dict["transition_group_id"]]
                      if unique_peptide_id == trgroup_id:
                          row_to_write = row
                          row_to_write += [selected_ids_dict[f_id].peptide.run.get_id(), f]
                          writer.writerow(row_to_write)

        self._write_trafo_files()

        if len(yaml_outfile) > 0:
            import yaml
            myYaml = {"RawData" : [], "PeakGroupData" : [ outfile ],
                      "ReferenceRun" : self.transformation_collection.getReferenceRunID() }
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


# Detect outliers in "good" groups
def detect_outliers(multipeptides, aligned_fdr_cutoff, outlier_threshold_seconds):
    class Outliers:
        def __init__(self): pass
    o = Outliers
    o.outliers_changed = 0
    o.outlier_pg = 0
    outliers = 0
    for m in multipeptides:
      out = m.detect_outliers()
      if len(out) == 0: continue
      rts = [float(p.get_selected_peakgroup().get_normalized_retentiontime()) for p in m.get_peptides() if p.get_selected_peakgroup() is not None]
      outlier_rts = [float(p.get_selected_peakgroup().get_normalized_retentiontime()) for p in m.get_peptides() if p.get_run_id() in out]
      mean_wo_outliers = numpy.mean([r for r in rts if r not in outlier_rts])
      for outlier_idx,outlier in enumerate(outlier_rts):
        if abs(outlier - mean_wo_outliers) > outlier_threshold_seconds:
            outliers += 1
            o.outlier_pg += 1
            thispep = [pep for pep in m.get_peptides() if pep.get_run_id() == out[outlier_idx]][0]
            newpg = thispep.find_closest_in_iRT(mean_wo_outliers)
            print "Bad", out, rts, "-> o", outlier_rts, "m", mean_wo_outliers, thispep.get_selected_peakgroup().get_normalized_retentiontime(), \
                " => exch", newpg.get_normalized_retentiontime(), newpg.get_fdr_score(), "/", thispep.get_selected_peakgroup().get_fdr_score()
            if( #abs(newpg.get_normalized_retentiontime() - best_rt_diff) < rt_diff_cutoff and
                newpg.get_fdr_score() < aligned_fdr_cutoff):
                  thispep.get_selected_peakgroup().unselect_this_peakgroup()
                  newpg.select_this_peakgroup()
                  o.outliers_changed += 1
                  print "change!"
        else: pass
    o.nr_outliers = outliers
    print "Changed %s outliers out of %s outliers (all pg %s)" % (o.outliers_changed, o.outlier_pg, 0)
    return o

def estimate_aligned_fdr_cutoff(options, this_exp, multipeptides, fdr_range):
    print "Try to find parameters for target fdr %0.2f %%" % (options.target_fdr * 100)
    for aligned_fdr_cutoff in fdr_range:
        # do the alignment and annotate chromatograms without identified features
        # then perform an outlier detection over multiple runs
        # unselect all
        for m in multipeptides:
            for p in m.get_peptides():
                p.unselect_all()
        # now align
        options.aligned_fdr_cutoff = aligned_fdr_cutoff
        alignment = align_features(multipeptides, options.rt_diff_cutoff, options.fdr_cutoff, options.aligned_fdr_cutoff, options.method)
        est_fdr = this_exp.estimate_real_fdr(multipeptides, options.min_frac_selected).est_real_fdr
        print "Estimated FDR: %0.4f %%" % (est_fdr * 100), "at position aligned fdr cutoff ", aligned_fdr_cutoff
        if est_fdr > options.target_fdr:
            # Unselect the peptides again ...
            for m in multipeptides:
                for p in m.get_peptides():
                    p.unselect_all()
            return aligned_fdr_cutoff

class DReadFilter(object):
    def __init__(self, cutoff):
        self.cutoff = cutoff
    def __call__(self, row, header):
        return float(row[ header["d_score" ] ]) > self.cutoff

class ParamEst(object):
    """
    Parameter estimation object

    In a first step the percentage of decoys of all peakgroups at the target
    fdr is computed (which is then taken as the "aim"). For this "aim" of decoy
    percentage, the class will try to estimate an fdr_cutoff such that the
    percentage of decoy precursors in the final reported result will correspond
    to the "aim". 

    If the parameter min_runs (at initialization) is higher than 1, only
    precursors that are identified in min_runs above the fdr_cutoff will be
    reported.

    >>> p = ParamEst()
    >>> decoy_frac = p.compute_decoy_frac(multipeptides, target_fdr)
    >>> fdr_cutoff_calculated = p.find_iterate_fdr(multipeptides, decoy_frac)
    """

    def __init__(self, min_runs=1, verbose=False):
        self.verbose = verbose
        self.min_runs = min_runs

    def find_iterate_fdr(self, multipeptides, decoy_frac, recursion=0):

        # Starting guess
        start = 0.05 / (10**recursion)
        end = 1.0 / (10**recursion)
        stepsize = start

        if self.verbose: print "Recurse", recursion
        if recursion > 10:
            raise Exception("Recursed too much in FDR iteration.")

        decoy_pcnt = decoy_frac*100
        val_005 = self._calc_precursor_fr(multipeptides, (start+stepsize)/100.0 )*100
        val_1 = self._calc_precursor_fr(multipeptides, end/100.0 )*100
        if self.verbose: print "Decoy pcnt aim:", decoy_pcnt
        if self.verbose: print decoy_pcnt, val_1, val_005

        # Check if our computed value lies between 0.05% and 1% FDR cutoff
        if decoy_pcnt < val_005:
            return self.find_iterate_fdr(multipeptides, decoy_frac, recursion=recursion+1)
        elif decoy_pcnt > val_1:
            if self.verbose: print "choose larger step from 0.5 on"
            start = 0.5
            end = 100.0
            stepsize = 0.5
            if recursion > 1:
                raise Exception("Decreased start / end but the values was too large? Should never happen.")
        else:
            # All is fine, we are within the limits
            pass

        fdrrange = numpy.arange(start, end + 2*stepsize, stepsize) # add 2 extra steps for edge cases
        return self._find_iterate_fdr(multipeptides, decoy_frac, fdrrange)

    def _find_iterate_fdr(self, multipeptides, decoy_frac, fdrrange):
        decoy_pcnt = decoy_frac*100
        if self.verbose: print "mScore_cutoff", "Calc-precursor-FDR"
        for fdr in fdrrange:
            calc_fdr = self._calc_precursor_fr(multipeptides, fdr/100.0 )*100
            if self.verbose: print fdr, calc_fdr
            if calc_fdr > decoy_pcnt:
                break
            prev_fdr = fdr
            prev_calc_fdr = calc_fdr

        # The last value is extremely close to the true one
        if abs(calc_fdr - decoy_pcnt) < 1e-6:
            return fdr/100.0

        # We have run through without stopping
        if abs(prev_fdr - fdr) < 1e-6:
            raise Exception("Parameter estimation did not reach a high enough value")

        # Linear interpolation
        res = prev_fdr + (fdr-prev_fdr) * (decoy_pcnt-prev_calc_fdr)/(calc_fdr-prev_calc_fdr)
        return res/100.0

    def _calc_precursor_fr(self, multipeptides, target_fdr):
        """ Calculate how many of the *precursors* are decoy for a given cutoff.
        """
        min_runs = self.min_runs
        ## min_runs = 1
        allpg_cnt = 0
        alldecoypg_cnt = 0
        for mpep in multipeptides:
            count = 0
            decoy = False
            for pep in mpep.get_peptides():
                if pep.get_best_peakgroup().get_fdr_score() < target_fdr:
                    count += 1
                if pep.get_decoy():
                    decoy = True
            if count >= min_runs:
                allpg_cnt += 1
            if decoy and count >= min_runs:
                alldecoypg_cnt += 1
        return alldecoypg_cnt *1.0 / allpg_cnt

    def compute_decoy_frac(self, multipeptides, target_fdr):
        """ Calculate how many of the *peakgroups* are decoy for a given cutoff.
        """
        allpg_cnt = 0
        alldecoypg_cnt = 0
        for mpep in multipeptides:
            for pep in mpep.get_peptides():
                if pep.get_best_peakgroup().get_fdr_score() < target_fdr:
                    allpg_cnt += 1
                    if pep.get_decoy():
                        alldecoypg_cnt += 1


        decoy_frac = alldecoypg_cnt *1.0 / allpg_cnt
        return decoy_frac

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program will select all peakgroups below the FDR cutoff in all files and try to align them to each other."
    usage += "\nIf only one file is given, it will act as peakgroup selector (best by m_score)" + \
            "\nand will apply the provided FDR cutoff."

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", required=True, nargs = '+', help = 'A list of mProphet output files containing all peakgroups (use quotes around the filenames)')
    parser.add_argument("--out", dest="outfile", required=True, default="feature_alignment_outfile", help="Output file with filtered peakgroups for quantification (only works for OpenSWATH)")
    parser.add_argument("--out_matrix", dest="matrix_outfile", default="", help="Matrix containing one peak group per row")
    parser.add_argument("--out_ids", dest="ids_outfile", default="", help="Id file only containing the ids")
    parser.add_argument("--out_meta", dest="yaml_outfile", default="", help="Outfile containing meta information, e.g. mapping of runs to original directories")
    parser.add_argument("--fdr_cutoff", dest="fdr_cutoff", default=0.01, type=float, help="FDR cutoff to use, default 0.01", metavar='0.01')
    parser.add_argument("--max_rt_diff", dest="rt_diff_cutoff", default=30, type=float, help="Maximal difference in RT for two aligned features", metavar='30')
    parser.add_argument("--max_fdr_quality", dest="aligned_fdr_cutoff", default=-1.0, help="Quality cutoff to still consider a feature for alignment (in FDR) - it is possible to give a range in the format lower,higher+stepsize,stepsize - e.g. 0,0.31,0.01 (-1 will set it to fdr_cutoff)", metavar='-1')
    parser.add_argument("--frac_selected", dest="min_frac_selected", default=0.0, type=float, help="Do not write peakgroup if selected in less than this fraction of runs (range 0 to 1)", metavar='0')
    parser.add_argument('--method', default='best_overall', help="Which method to use for the clustering (best_overall, best_cluster_score or global_best_cluster_score, global_best_overall). The global option will also move peaks which are below the selected FDR threshold.")
    parser.add_argument('--file_format', default='openswath', help="Which input file format is used (openswath or peakview)")

    experimental_parser = parser.add_argument_group('experimental options')

    experimental_parser.add_argument('--use_dscore_filter', action='store_true', default=False)
    experimental_parser.add_argument("--dscore_cutoff", default=1.96, type=float, help="Quality cutoff to still consider a feature for alignment using the d_score: everything below this d-score is discarded", metavar='1.96')
    experimental_parser.add_argument("--nr_high_conf_exp", default=1, type=int, help="Number of experiments in which the peptide needs to be identified with high confidence (e.g. above fdr_curoff)", metavar='1')
    experimental_parser.add_argument("--readmethod", dest="readmethod", default="minimal", help="Read full or minimal transition groups (minimal,full)")
    experimental_parser.add_argument("--outlier_thresh", dest="outlier_threshold_seconds", default=30, type=float, help="Everything below this threshold (in seconds), a peak will not be considered an outlier", metavar='30')
    experimental_parser.add_argument('--remove_outliers', action='store_true', default=False)
    experimental_parser.add_argument('--realign_runs', action='store_true', default=False, help="Tries to re-align runs based on their true RT (instead of using the less accurate iRT values by computing a spline against a reference run)")
    experimental_parser.add_argument('--use_scikit', action='store_true', default=False, help="Use datasmooth from scikit instead of R to re-align runs (needs to be installed)")
    experimental_parser.add_argument('--use_linear', action='store_true', default=False, help="Use linear run alignment")
    experimental_parser.add_argument('--use_external_r', action='store_true', default=False, help="Use external R call for alignment (instead of rpy2)")
    experimental_parser.add_argument("--tmpdir", dest="tmpdir", default="/tmp/", help="Temporary directory")
    experimental_parser.add_argument("--alignment_score", dest="alignment_score", default=0.0001, type=float, help="Minimal score needed for a feature to be considered for alignment between runs", metavar='0.0001')
    experimental_parser.add_argument("--target_fdr", dest="target_fdr", default=-1, type=float, help="If parameter estimation is used, which target FDR should be optimized for. If set to lower than 0, parameter estimation is turned off.", metavar='0.01')

    args = parser.parse_args(sys.argv[1:])

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
    import time

    readfilter = ReadFilter()
    if options.use_dscore_filter:
        readfilter = DReadFilter(float(options.dscore_cutoff))

    # Read the files
    start = time.time()
    reader = SWATHScoringReader.newReader(options.infiles, options.file_format, options.readmethod, readfilter)
    runs = reader.parse_files(options.realign_runs)
    # Create experiment
    this_exp = Experiment()
    this_exp.set_runs(runs)
    print("Reading the input files took %ss" % (time.time() - start) )

    # Map the precursors across multiple runs, determine the number of
    # precursors in all runs without alignment.
    start = time.time()
    multipeptides = this_exp.get_all_multipeptides(options.fdr_cutoff, verbose=True)
    print("Mapping the precursors took %ss" % (time.time() - start) )

    if options.target_fdr > 0:
        ### Do parameter estimation
        start = time.time()
        print "-"*35
        print "Do Parameter estimation"
        p = ParamEst(min_runs=options.nr_high_conf_exp,verbose=True)
        decoy_frac = p.compute_decoy_frac(multipeptides, options.target_fdr)
        print "Found target decoy fraction overall %0.4f%%" % (decoy_frac*100)
        try:
            fdr_cutoff_calculated = p.find_iterate_fdr(multipeptides, decoy_frac)
        except UnboundLocalError:
            raise Exception("Could not estimate FDR accurately!")

        if fdr_cutoff_calculated > options.target_fdr:
            #### Re-read the multipeptides with the new cutoff ... !
            multipeptides = this_exp.get_all_multipeptides(fdr_cutoff_calculated, verbose=True)
            print "Re-parse the files!"
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
        print "Using an FDR cutoff of %0.4f%%" % (fdr_cutoff_calculated*100)
        print "For the aligned values, use a cutoff of %0.4f%%" % (options.aligned_fdr_cutoff)
        print("Parameter estimation took %ss" % (time.time() - start) )
        print "-"*35

    # If we want to align runs
    if options.realign_runs:
        start = time.time()
        spl_aligner = SplineAligner(options.alignment_score,
                                   options.use_scikit,
                                   options.use_linear,
                                   options.use_external_r, options.tmpdir)
        tcoll = spl_aligner.rt_align_all_runs(this_exp, multipeptides)
        this_exp.transformation_collection = tcoll

        print("Aligning the runs took %ss" % (time.time() - start) )

    try:
        options.aligned_fdr_cutoff = float(options.aligned_fdr_cutoff)
    except ValueError:
        # We have a range, since we trust the input we dont parse it very much ...
        exec("fdr_range = numpy.arange(%s)" % options.aligned_fdr_cutoff)
        options.aligned_fdr_cutoff = estimate_aligned_fdr_cutoff(options, this_exp, multipeptides, fdr_range)

    print "Will calculate with aligned_fdr cutoff of", options.aligned_fdr_cutoff
    start = time.time()
    alignment = AlignmentAlgorithm().align_features(multipeptides, options.rt_diff_cutoff, options.fdr_cutoff, options.aligned_fdr_cutoff, options.method)

    print("Re-aligning peak groups took %ss" % (time.time() - start) )
    if options.remove_outliers:
      outlier_detection = detect_outliers(multipeptides, options.aligned_fdr_cutoff, options.outlier_threshold_seconds)
    else: outlier_detection = None

    # Filter by high confidence (e.g. keep only those where enough high confidence IDs are present)
    for mpep in multipeptides:
        # check if we have found enough peakgroups which are below the cutoff
        count = 0
        for pg in mpep.get_selected_peakgroups():
            if pg.get_fdr_score() < options.fdr_cutoff:
                count += 1
        if count < options.nr_high_conf_exp:
            for p in mpep.get_peptides():
                p.unselect_all()

    # print statistics, write output
    start = time.time()
    this_exp.print_stats(multipeptides, alignment, outlier_detection, options.fdr_cutoff, options.min_frac_selected, options.nr_high_conf_exp)
    this_exp.write_to_file(multipeptides, options)
    print("Writing output took %ss" % (time.time() - start) )

if __name__=="__main__":
    options = handle_args()
    main(options)

