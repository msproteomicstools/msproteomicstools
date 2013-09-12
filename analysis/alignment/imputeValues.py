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
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import AlignmentExperiment, Multipeptide
from sys import stdout

from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import AlignmentExperiment as Experiment 

class ImputeValuesHelper(object):
    """
    Static object with some helper methods.
    """

    @staticmethod
    def select_correct_swath(swath_chromatograms, mz):
        swath_window_low = int(mz / 25) *25
        swath_window_high = int(mz / 25) *25 + 25
        res = {}
        for k,v in swath_chromatograms.iteritems():
            selected = [vv for prec_mz,vv in v.iteritems() if prec_mz > swath_window_low and prec_mz < swath_window_high]
            if len(selected) == 1: 
                res[k] = selected[0]
        return res

    @staticmethod
    def convert_to_this(orig_runid, target_runid, ref_id, rt, transformation_collection_):
        """ Convert a retention time into one of the target RT space.
        
        Using the transformation collection
        """
        normalized_space_rt = transformation_collection_.getTransformation(orig_runid, ref_id).predict( [rt] )[0]
        return transformation_collection_.getTransformation(ref_id, target_runid).predict( [normalized_space_rt] )[0]

def read_chromatograms(options, peakgroups_file, trafo_fnames):

    # options.outfile = "/tmp/testout.csv"
    # trafo_fnames = ['analysis/alignment/strep_align/Strep0_Repl1_R02/transformation-0_0-0_1.tr', 'analysis/alignment/strep_align/Strep0_Repl2_R02/transformation-0_1-0_1.tr', 'analysis/alignment/strep_align/Strep10_Repl1_R02/transformation-0_2-0_1.tr', 'analysis/alignment/strep_align/Strep10_Repl2_R02/transformation-0_3-0_1.tr']

    # trafo_fnames = ["strep_align/Strep0_Repl1_R02/transformation-0_0-0_0.tr", "strep_align/Strep0_Repl2_R02/transformation-0_1-0_0.tr", "strep_align/Strep10_Repl1_R02/transformation-0_2-0_0.tr", "strep_align/Strep10_Repl2_R02/transformation-0_3-0_0.tr" ]
    # options.outfile = "../../output_010_Repl12_R02/out_reduced.txt"
    # peakgroups_file = "../../output_010_Repl12_R02/out_reduced.txt"
    # options.file_format = "openswath"
    # options.fdr_cutoff = "0.01"
    # options.realign_runs = False

    fdr_cutoff_all_pg = 1.0
    reader = SWATHScoringReader.newReader([peakgroups_file], options.file_format, readmethod="complete")
    new_exp = Experiment()
    new_exp.runs = reader.parse_files()
    multipeptides = new_exp.get_all_multipeptides(fdr_cutoff_all_pg, verbose=False)

    transformation_collection_ = TransformationCollection()
    for filename in trafo_fnames:
      # current_id = current_run.get_id()
      # filename = os.path.join(os.path.dirname(current_run.orig_filename), "transformation-%s-%s.tr" % (current_id, ref_id) )
      transformation_collection_.readTransformationData(filename)

    # Read the datapoints and perform the smoothing
    transformation_collection_.initialize_from_data(reverse=True)

    # Read the mzML files and store them
    # TODO abstract this away to an object
    swath_chromatograms = {}
    for filename in trafo_fnames:
        # get the run id
        f = open(filename, "r")
        header = f.next().split("\t")
        f.close()
        all_swathes = {}
        runid = header[1]
        import glob
        dname = os.path.dirname(filename)
        files = glob.glob(os.path.join(dname + "/*.mzML") )
        for f in files:
            import pymzml
            run = pymzml.run.Reader(f, build_index_from_scratch=True)
            first = run.next()
            mz = first['precursors'][0]['mz']
            all_swathes[ int(mz) ] = run
        swath_chromatograms[ runid ] = all_swathes

    if options.dry_run:
        print "Dry Run only"
        print "Found multipeptides:", len(multipeptides)
        print "Found swath chromatograms:", len(swath_chromatograms), [v.keys() for k,v in swath_chromatograms.iteritems()]
        return [], []

    return analyze_multipeptides(new_exp, multipeptides, swath_chromatograms, transformation_collection_)
    
def analyze_multipeptides(new_exp, multipeptides, swath_chromatograms, transformation_collection_):
    # Go through all aligned peptides
    run_ids = [r.get_id() for r in new_exp.runs]
    imputations = 0
    imputation_succ = 0
    peakgroups = 0
    for m in multipeptides:
        line = [m.get_id()]
        selected_pg = [p.peakgroups[0] for p in m.get_peptides() if len(p.peakgroups)==1 ] 
        current_mz = float(selected_pg[0].get_value("m.z"))
        # return [p.get_selected_peakgroup() for p in self.peptides.values() if p.get_selected_peakgroup() is not None]
        # print "selected pgs ... ", len(selected_pg)
        for rid in run_ids:
            pg = None
            peakgroups += 1
            #print m
            #print m.get_peptides(), rid, rid in m.get_peptides()
            if m.has_peptide(rid):
                if len( m.get_peptide(rid).peakgroups) > 1:
                    print "found more than one pg"
                    print m._peptides[rid].peakgroups
                    for pg in  m.peptides[rid].peakgroups:
                        print pg.row
                assert len( m.get_peptide(rid).peakgroups) == 1
                pg = m.get_peptide(rid).peakgroups[0]
                # Also select this pg
                pg.select_this_peakgroup()
            if pg is None:
                # fill up the hole if we have an NA here!
                imputations += 1

                current_run = [r for r in new_exp.runs if r.get_id() == rid][0]

                border_l, border_r = determine_integration_border(new_exp, selected_pg, 
                    swath_chromatograms, rid, transformation_collection_)
                res = integrate_chromatogram(selected_pg[0], current_run, swath_chromatograms, current_mz,
                                             border_l, border_r)
                if res != "NA": 
                    imputation_succ += 1
                    # print "result is", res.get_value("Intensity"), numpy.log10(res.get_value("Intensity"))
                    p = GeneralPrecursor(res.get_value("transition_group_id"), rid)
                    # Also select this pg
                    res.select_this_peakgroup()
                    p.add_peakgroup(res)
                    m.add_peptide(rid, p)
    print "Imputations:", imputations, "Successful:", imputation_succ
    print "Peakgroups:", peakgroups
    return new_exp, multipeptides 

def determine_integration_border(new_exp, selected_pg, swath_chromatograms, rid, transformation_collection_):
    """Determine the optimal integration border by taking the mean of all other peakgroup boundaries.

    Args:
        new_exp(AlignmentExperiment): experiment containing the multipeptides
        selected_pg(list(GeneralPeakGroup)): list of selected peakgroups (e.g. those passing the quality threshold)
        swath_chromatograms(dict): containing the objects pointing to the original chrom mzML
        rid(String): current run id
        transformation_collection_(.TransformationCollection): specifying how to transform between retention times of different runs

    Returns:
        A tuple of (left_integration_border, right_integration_border) in the retention time space of the _reference_ run
    """
    current_run = [r for r in new_exp.runs if r.get_id() == rid][0]
    ref_id = transformation_collection_.reference_run_id

    pg_lefts = []
    pg_rights = []
    for pg in selected_pg:

        rwidth = float(pg.get_value("rightWidth"))
        lwidth = float(pg.get_value("leftWidth"))
        this_run_rwidth = ImputeValuesHelper.convert_to_this(pg.peptide.run.get_id(),
            current_run.get_id(), ref_id, rwidth, transformation_collection_)
        this_run_lwidth = ImputeValuesHelper.convert_to_this(pg.peptide.run.get_id(),
            current_run.get_id(), ref_id, lwidth, transformation_collection_)

        pg_lefts.append(this_run_lwidth)
        pg_rights.append(this_run_rwidth)

        # print pg.peptide.run.get_id(), "look at peakgroup %0.4f" % this_run_lwidth,\
        #       "\t",this_run_rwidth , "\t",this_run_rwidth-this_run_lwidth, "\tscore:", \
        #       -numpy.log10(float(pg.get_value("m_score"))),  "intensity:", \
        #       numpy.log10(float(pg.get_value("Intensity"))), pg.get_value("Sequence")

    # print "overall mean, std", numpy.mean(pg_lefts), numpy.std(pg_lefts), numpy.mean(pg_rights), numpy.std(pg_rights)

    integration_right = numpy.mean(pg_rights)
    integration_left = numpy.mean(pg_lefts)

    std_warning_level = 20
    if numpy.std(pg_rights) > std_warning_level or numpy.std(pg_lefts) > std_warning_level: 
        pass
        # TODO : what if they are not consistent ?
        # print "Std is too large", pg_rights, pg_lefts
        # print "overall mean, std", numpy.mean(pg_lefts), numpy.std(pg_lefts), numpy.mean(pg_rights), numpy.std(pg_rights)

    return integration_left, integration_right

def integrate_chromatogram(template_pg, current_run, swath_chromatograms, current_mz, 
                           left_start, right_end):
    """ Integrate a chromatogram from left_start to right_end and store the sum.

    Args:
        pg(GeneralPeakGroup): A template peakgroup from which to construct the new peakgroup
        current_run(SWATHScoringReader.Run): current run where the missing value occured
        current_mz(float): m/z value of the current precursor
        left_start(float): retention time for integration (left border)
        right_end(float): retention time for integration (right border)

    Returns:
        A new GeneralPeakGroup which contains the new, integrated intensity for this run (or "NA" if no chromatograms could be found).

    Create a new peakgroup from the old pg and then store the integrated intensity.
    """

    # get the chromatograms
    current_rid  = current_run.get_id()
    correct_swath = ImputeValuesHelper.select_correct_swath(swath_chromatograms, current_mz)
    chrom_ids = template_pg.get_value("aggr_Fragment_Annotation").split(";")
    if not correct_swath.has_key(current_rid):
        return "NA"

    # Create new peakgroup by copying the old one
    newrow = ["NA" for ele in template_pg.row]
    newpg = GeneralPeakGroup(newrow, current_run, template_pg.peptide)
    for element in ["transition_group_id", "decoy", "Sequence", "FullPeptideName", "Charge", "ProteinName", "nr_peaks"]:
        newpg.set_value(element, template_pg.get_value(element))
    newpg.set_value("transition_group_record", template_pg.get_value("transition_group_id") + "_%s" % current_rid)
    newpg.set_value("run_id", current_rid)
    newpg.set_value("RT", (left_start + right_end) / 2.0 )
    newpg.set_value("leftWidth", left_start)
    newpg.set_value("rightWidth", right_end)

    integrated_sum = 0
    allchroms = correct_swath[current_rid]
    for chrom_id in chrom_ids:
        chromatogram = allchroms[chrom_id]
        if chromatogram is None:
            print "chromatogram is None"
            continue
        integrated_sum += sum( [p[1] for p in chromatogram.peaks if p[0] > left_start and p[0] < right_end ])
        # print integrated_sum, "chromatogram", sum(p[1] for p in chromatogram.peaks), \
        # "integrated from \t%s\t%s in run %s" %( left_start, right_end, current_rid)

    newpg.set_value("Intensity", integrated_sum)
    return newpg

def write_out(new_exp, multipeptides, outfile):
    # write out the complete original files 
    writer = csv.writer(open(outfile, "w"), delimiter="\t")
    header_first = new_exp.runs[0].header
    for run in new_exp.runs:
        assert header_first == run.header
    writer.writerow(header_first)

    print "len multipeptides"
    for m in multipeptides:
        # selected_peakgroups = [p.peakgroups[0] for p in m.get_peptides()]
        # if (len(selected_peakgroups)*2.0 / len(new_exp.runs) < fraction_needed_selected) : continue
        for p in m.get_peptides():
            selected_pg = p.peakgroups[0]
            if selected_pg is None: 
                continue
            row_to_write = selected_pg.row
            row_to_write += [selected_pg.run.get_id(), selected_pg.run.orig_filename]
            writer.writerow(row_to_write)

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program will impute missing values"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", nargs = '+', required=True, help = 'A list of mProphet output files containing all peakgroups (use quotes around the filenames)')
    parser.add_argument("--peakgroups_infile", dest="peakgroups_infile", required=True, help="Infile containing peakgroups (outfile from feature_alignment.py)")
    parser.add_argument("--out", dest="output", required=True, help="Output file with imputed values")
    parser.add_argument('--file_format', default='openswath', help="Which input file format is used (openswath or peakview)")
    parser.add_argument('--dry_run', action='store_true', default=False, help="Perform a dry run only")

    experimental_parser = parser.add_argument_group('experimental options')

    args = parser.parse_args(sys.argv[1:])
    return args

def main(options):
    import time
    new_exp, multipeptides = read_chromatograms(options, options.peakgroups_infile, options.infiles)
    if options.dry_run: return
    write_out(new_exp, multipeptides, options.output)

if __name__=="__main__":
    options = handle_args()
    main(options)

