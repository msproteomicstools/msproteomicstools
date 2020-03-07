#!/usr/bin/env python
# -*- coding: utf-8  -*-
import os, sys, csv, time

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
from msproteomicstoolslib.algorithms.alignment.FDRParameterEstimation import ParamEst
from msproteomicstoolslib.algorithms.PADS.MinimumSpanningTree import MinimumSpanningTree

import numpy as np
from analysis.alignment.feature_alignment import Experiment
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.format.SWATHScoringMapper import MSfileRunMapping, getPrecursorTransitionMapping
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
import matplotlib.pyplot as plt
fig = plt.figure()
plt.close()
from pyopenms import OnDiscMSExperiment
import msproteomicstoolslib.math.Smoothing as smoothing

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


def extractXIC_group(mz, chromIndices):
    """ Extract chromatograms for chromatrogram indices """
    XIC_group = [None]*len(chromIndices)
    for i in range(len(chromIndices)):
        # Chromatogram is a tuple (time_array, intensity_array) of numpy array
        XIC_group[i] = mz.getChromatogram(chromIndices[i]).get_peaks()
    return XIC_group

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
precursor_to_transitionID = getPrecursorTransitionMapping(infiles[0])

# Establish connection to mzML files
# TODO Add the mz accessor to Run object
chrom_file_accessor = {}
for run in runs:
    mz = OnDiscMSExperiment()
    mzml_file = MStoFeature[run.get_openswath_filename()][0]
    mz.openFile(mzml_file)
    chrom_file_accessor[run.get_id()] = mz

# Get precursor to chromatogram indices
run_chromIndex_map = {}
for run in runs:
    meta_data = chrom_file_accessor[run.get_id()].getMetaData()
    precursor_to_chromIndex = {}
    # Initialize chromatogram indices with -1
    for prec in precursor_to_transitionID:
        precursor_to_chromIndex[prec] = [-1]*len(precursor_to_transitionID[prec])
    # Iterate through Native IDs in mzML file
    for i in range(meta_data.getNrChromatograms()):
        nativeID = int(meta_data.getChromatogram(i).getNativeID())
        # Check if Native ID matches to any of our transition IDs.
        for prec in precursor_to_transitionID:
            if nativeID in precursor_to_transitionID[prec]:
                # Find the index of transition ID and insert chromatogram index
                index = precursor_to_transitionID[prec].index(nativeID)
                precursor_to_chromIndex[prec][index] = i
    run_chromIndex_map[run.get_id()] = precursor_to_chromIndex

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

# Pairwise global alignment
spl_aligner = SplineAligner(alignment_fdr_threshold = 0.05, smoother="lowess", experiment=this_exp)
spl_aligner.initialize_transformation_error()
# Initialize XIC smoothing function
sm = smoothing.getXIC_SmoothingObj(smoother = "sgolay", kernelLen = 11, polyOrd = 4)
def smoothXICs(XIC_Group, sm):
    XIC_Group_sm = []
    for XIC in XIC_Group:
        sm.initialize(XIC[0], XIC[1])
        XIC_sm = sm.smooth(XIC[0], XIC[1])
        XIC_Group_sm.append(XIC_sm)
    return XIC_Group_sm

# Use multiprocessing tool for extracting chromatograms
prec_id = 9719 #, 9720
RSEdistFactor = 4
for prec_id in precursor_to_transitionID:
    refrun = reference_run[prec_id]
    eXps = list(set(runs) - set([refrun]))
    # Extract XICs from reference run
    refrun_id = refrun.get_id()
    chrom_ids = run_chromIndex_map[refrun_id][prec_id]
    XICs_ref = extractXIC_group(chrom_file_accessor[refrun_id], chrom_ids)
    XICs_ref_sm = smoothXICs(XICs_ref, sm)
    # Iterate through all other runs and align them to the reference run
    for eXprun in eXps:
        ## Extract XICs from experiment run
        eXprun_id = eXprun.get_id()
        chrom_ids = run_chromIndex_map[eXprun_id][prec_id]
        XICs_eXp = extractXIC_group(chrom_file_accessor[eXprun_id], chrom_ids)
        XICs_eXp_sm = smoothXICs(XICs_eXp, sm)
        # Get time component
        t_ref = XICs_ref_sm[0][0]
        t_eXp = XICs_eXp_sm[0][0]
        ## Set up constraints for penalizing similarity matrix
        # Get constraining end-points using Global alignment
        globalAligner = spl_aligner.rt_align_pair(refrun, eXprun, multipeptides)
        B1p = globalAligner.predict([ t_ref[0] ])
        B2p = globalAligner.predict([ t_ref[-1] ])
        # Get width of the constraining window
        adaptiveRT = RSEdistFactor * spl_aligner.transformation_error.transformations[refrun_id][eXprun_id][0]
        samplingTime = (t_ref[-1] - t_ref[0]) / (len(t_ref) - 1)
        noBeef = np.ceil(adaptiveRT/samplingTime)
        ## Get intensity values to build similarity matrix
        intensityList_ref = np.array([xic[1] for xic in XICs_ref_sm])
        intensityList_eXp = np.array([xic[1] for xic in XICs_eXp_sm])
        intensityList_ref = np.ascontiguousarray(intensityList_ref)
        intensityList_eXp = np.ascontiguousarray(intensityList_eXp)

    
    
# Iterate through each precursor and get the alignment.
    #Get smoothing from SplineAligner. Check if smoothing is already present.
    #calculate AdaptiveRT
    #getAlignObj function that calls DIAlignPy
    # SplineAligner.py _spline_align_runs(self, bestrun, run, multipeptides)

fig = plt.figure()
XIC_Group = XICs_ref_sm
for k in range(len(XIC_Group)):
    x = XIC_Group[k][0]
    y = XIC_Group[k][1]
    plt.plot(x, y)

plt.ylabel('Intensity')
plt.show()
plt.close()