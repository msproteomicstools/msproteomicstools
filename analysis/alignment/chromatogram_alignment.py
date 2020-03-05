#!/usr/bin/env python
# -*- coding: utf-8  -*-
import os, sys, csv, time
import numpy
import argparse
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.version import __version__ as version
from msproteomicstoolslib.math.chauvenet import chauvenet
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection, LightTransformationData
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment
from msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm import AlignmentAlgorithm
from msproteomicstoolslib.algorithms.alignment.AlignmentMST import getDistanceMatrix, TreeConsensusAlignment
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import write_out_matrix_file, addDataToTrafo
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
from msproteomicstoolslib.algorithms.alignment.FDRParameterEstimation import ParamEst
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree

from analysis.alignment.feature_alignment import Experiment
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.format.SWATHScoringMapper import MSfileRunMapping, getPrecursorTransitionMapping
#from analysis.alignment.chromatogram_alignment import get_reference_run, get_precursor_reference_run, get_multipeptide_reference_run

# Get a single reference run
def get_reference_run(experiment, multipeptides, alignment_fdr_threshold = 0.05):
    """
    Returns a dectionary with precursor id as key and run as value.
    Single reference run for each precursor.
    """
    best_run = experiment.determine_best_run(alignment_fdr_threshold = 0.05)
    reference_run = {}
    precursor_ids = set()
    for i in range(len(multipeptides)):
        # Get precursor groups from each multipeptide
        prec_groups =  multipeptides[i].getPrecursorGroups()
        for prec_group in prec_groups:
            # Get precursors from a precursor group
            precs = prec_group.getAllPrecursors()
            for prec in precs:
                # Get all precursor ids
                precursor_ids.add(prec.get_id())
    # Get a sorted list
    precursor_ids = sorted(precursor_ids)
    # Assign best_run as reference run for each precursor_id
    reference_run = dict.fromkeys(precursor_ids , best_run)
    return reference_run    

# Get reference run for each precursor
def get_precursor_reference_run(multipeptides, alignment_fdr_threshold = 0.05):
    """
    Returns a dectionary with precursor id as key and run as value.
    Precursors may have different reference runs based on FDR score of associated best peak group.
    """
    reference_run = {}
    for i in range(len(multipeptides)):
        # Get precursor groups from each multipeptide
        prec_groups =  multipeptides[i].getPrecursorGroups()
        for prec_group in prec_groups:
            # Get precursors from a precursor group
            precs = prec_group.getAllPrecursors()
            max_fdr = 1.0
            for prec in precs:
                # Get all precursor ids
                prec_id = prec.get_id()
                # Get best FDR value for the precursor
                cur_fdr = prec.get_best_peakgroup().get_fdr_score()
                if cur_fdr <= alignment_fdr_threshold and cur_fdr < max_fdr:
                    max_fdr = cur_fdr
                    # Make Run of the current precursor as the reference run
                    reference_run[prec_id] = prec.getRun()
            if prec_id not in reference_run:
                # No peak group has FDR lower than alignment_fdr_threshold
                reference_run[prec_id] = None
    return reference_run

# Get reference run for each precursor-group
def get_multipeptide_reference_run(multipeptides, alignment_fdr_threshold = 0.05):
    """
    Returns a dectionary with precursor id as key and run as value.
    Precursors may have different reference runs based on FDR score of associated multipeptide's best peak group.
    """
    reference_run = {}
    for i in range(len(multipeptides)):
        # Get precursor groups from each multipeptide
        prec_groups =  multipeptides[i].getPrecursorGroups()
        max_fdr = 1.0
        refRun = None
        for prec_group in prec_groups:
            precs = prec_group.getAllPrecursors()
            cur_fdr = prec_group.getOverallBestPeakgroup().get_fdr_score()
            if cur_fdr <= alignment_fdr_threshold and cur_fdr < max_fdr:
                max_fdr = cur_fdr
                refRun = prec_group.run_
            for prec in precs:
                # Get all precursor ids
                prec_id = prec.get_id()
                # Make Run of the current precursor as the reference run
                reference_run[prec_id] = refRun
    return reference_run


infiles = ["../DIAlignR/inst/extdata/osw/merged.osw"]
chromatograms = ["../DIAlignR/inst/extdata/mzml/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML",
 "../DIAlignR/inst/extdata/mzml/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML",
 "../DIAlignR/inst/extdata/mzml/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"]

readfilter = ReadFilter()
file_format = 'openswath'
readmethod = "minimal"

reader = SWATHScoringReader.newReader(infiles, file_format, readmethod, readfilter,
                                        enable_isotopic_grouping = False, read_cluster_id = False)
# reader.map_infiles_chromfiles(chromatograms)
runs = reader.parse_files()
MStoFeature = MSfileRunMapping(chromatograms, runs)
# MStoFeature[runs[0].get_openswath_filename()][0]
precursor_mapping = getPrecursorTransitionMapping(infiles[0])

this_exp = Experiment()
this_exp.set_runs(runs)
start = time.time()
fdr_cutoff = 0.05
multipeptides = this_exp.get_all_multipeptides(fdr_cutoff, verbose=False, verbosity= 10)
print("Mapping the precursors took %0.2fs" % (time.time() - start) )

# Reference based alignment
reference_run = get_reference_run(this_exp, multipeptides, alignment_fdr_threshold = 0.05)
reference_run = get_precursor_reference_run(multipeptides, alignment_fdr_threshold = 0.05)
reference_run = get_multipeptide_reference_run(multipeptides, alignment_fdr_threshold = 0.05)

# Iterate through each precursor and get the alignment.
    #Get precursor id
    #Get associated transition id
    #Extract XICs using pyopenms
    #smooth XICs using pyopenms
    #Get smoothing from SplineAligner. Check if smoothing is already present.
    #calculate AdaptiveRT
    #getAlignObj function that calls DIAlignPy
    # SplineAligner.py _spline_align_runs(self, bestrun, run, multipeptides)


for i in range(len(multipeptides)):
    for run in runs:
        id = run.get_id()
        prec = multipeptides[0].getPrecursorGroup(id).getAllPrecursors()
        



options = handle_args()

# python ./analysis/alignment/chromatogram_alignment.py --in file1_input.csv file2_input.csv file3_input.csv --out aligned.csv  --method best_overall --max_rt_diff 90 --target_fdr 0.01 --max_fdr_quality 0.05

# multipeptides[i].get_peptides()[j].get_run()
#                                   .get_precursors()

getAllPeptides # Returns a list of Precursor from each run. How come there are only three?
getPrecursorGroup
getPrecursorGroups


runs[0].get_id() # 125704171604355508
runs[1].get_id() # 6752973645981403097
runs[2].get_id() # 2234664662238281994

multipeptides[0].getPrecursorGroup(6752973645981403097)
multipeptides[0].getPrecursorGroup(2234664662238281994).getOverallBestPeakgroup().get_fdr_score()
multipeptides[0].getPrecursorGroup(2234664662238281994).getAllPrecursors()[0].get_best_peakgroup().get_fdr_score()

id = multipeptides[0].getPrecursorGroup(2234664662238281994).getAllPrecursors()[0].get_id()
id = 9719
id = 9720
chromid = precursor_mapping[id]

run = pymzml.run.Reader(chromatograms[0], build_index_from_scratch=True)
run.info["seekable"]
run.info["offsets"].items()
run[str(chromid[0])].peaks

getOverallBestPeakgroup
getAllPrecursors
getAllPeakgroups