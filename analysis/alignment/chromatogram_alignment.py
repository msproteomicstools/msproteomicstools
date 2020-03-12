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
import pandas as pd
from analysis.alignment.feature_alignment import Experiment
from msproteomicstoolslib.format.SWATHScoringReader import *
from msproteomicstoolslib.format.SWATHScoringMapper import MSfileRunMapping, getPrecursorTransitionMapping
from msproteomicstoolslib.algorithms.alignment.SplineAligner import SplineAligner
import matplotlib.pyplot as plt
fig = plt.figure()
plt.close()
from pyopenms import OnDiscMSExperiment
import msproteomicstoolslib.math.Smoothing as smoothing
from DIAlignPy import AffineAlignObj, alignChromatogramsCpp
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection, LightTransformationData
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import addDataToTrafo

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


def initialize_transformation():
    # Initialize objects to get transformation
    try:
        from msproteomicstoolslib.cython.LightTransformationData import CyLightTransformationData
        if optimized_cython:
            tr_data = CyLightTransformationData()
        else:
            tr_data = LightTransformationData()
    except ImportError:
        print("WARNING: cannot import CyLightTransformationData, will use Python version (slower).")
        tr_data = LightTransformationData()
    return tr_data

def get_pair_trafo(tr_data, spl_aligner, run_0, run_1, multipeptides,
                    max_rt_diff = 30, force = False, optimized_cython = False):
    """
    Returns the smoothing object that aligns run1 against run0.
    """
    run0_id = run_0.get_id()
    run1_id = run_1.get_id()
    # Add pairwise transformation if not found
    tmp = tr_data.trafo.get(run1_id, {})
    if not tmp.get(run0_id):
        addDataToTrafo(tr_data, run_0, run_1, spl_aligner, multipeptides,
                        realign_method = spl_aligner.smoother, max_rt_diff = max_rt_diff, force=force)
    return tr_data.getTrafo(run1_id, run0_id), tr_data.getStdev(run1_id, run0_id)

def extractXIC_group(mz, chromIndices):
    """ Extract chromatograms for chromatrogram indices """
    XIC_group = [None]*len(chromIndices)
    for i in range(len(chromIndices)):
        # Chromatogram is a tuple (time_array, intensity_array) of numpy array
        # mz.getChromatogramById(chromIndices[i]).getIntensityArray()
        # mz.getChromatogramById(chromIndices[i]).getTimeArray()
        XIC_group[i] = mz.getChromatogram(chromIndices[i]).get_peaks()
    return XIC_group

def get_aligned_time(t, indices, skipValue = -1):
    t_aligned = [None]*len(indices)
    for i in range(len(indices)):
        index = indices[i]
        if index != skipValue:
            t_aligned[i] = t[index]
    return t_aligned

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
invalid_chromIndex = -1
run_chromIndex_map = {}
for run in runs:
    meta_data = chrom_file_accessor[run.get_id()].getMetaData()
    precursor_to_chromIndex = {}
    # Initialize chromatogram indices with -1
    for prec in precursor_to_transitionID:
        precursor_to_chromIndex[prec] = [invalid_chromIndex]*len(precursor_to_transitionID[prec][1])
    # Iterate through Native IDs in mzML file
    for i in range(meta_data.getNrChromatograms()):
        nativeID = int(meta_data.getChromatogram(i).getNativeID())
        # Check if Native ID matches to any of our transition IDs.
        for prec in precursor_to_transitionID:
            if nativeID in precursor_to_transitionID[prec][1]:
                # Find the index of transition ID and insert chromatogram index
                index = precursor_to_transitionID[prec][1].index(nativeID)
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
tr_data = initialize_transformation()
# Initialize XIC smoothing function
sm = smoothing.getXIC_SmoothingObj(smoother = "sgolay", kernelLen = 11, polyOrd = 4)
def smoothXICs(XIC_Group, sm):
    XIC_Group_sm = []
    for XIC in XIC_Group:
        sm.initialize(XIC[0], XIC[1])
        XIC_sm = sm.smooth(XIC[0], XIC[1])
        XIC_Group_sm.append(XIC_sm)
    return XIC_Group_sm

# TODO Use multiprocessing tool for extracting chromatograms
RSEdistFactor = 4
# Calculate the aligned retention time for each precursor across all runs
prec_ids = list(precursor_to_transitionID.keys())
retention_time = {}
for run in runs:
    retention_time[run.get_id()] = [None]*len(prec_ids)

for i in range(len(prec_ids)):
    prec_id = prec_ids[i]
    refrun = reference_run.get(prec_id)
    if not refrun:
        print("The precursor {} doesn't have any associated reference run. Skipping!".format(prec_id))
        continue
    eXps = list(set(runs) - set([refrun]))
    # Extract XICs from reference run
    refrun_id = refrun.get_id()
    chrom_ids = run_chromIndex_map[refrun_id][prec_id]
    if invalid_chromIndex in chrom_ids:
        print("Can't extract XICs for precursor {} from run {}. Skipping!".format(prec_id, refrun_id))
        continue
    XICs_ref = extractXIC_group(chrom_file_accessor[refrun_id], chrom_ids)
    XICs_ref_sm = smoothXICs(XICs_ref, sm)
    # For each precursor, we need peptide_group_label and trgr_id
    peptide_group_label = precursor_to_transitionID[prec_id][0]
    refRT = refrun.getPrecursor(peptide_group_label, prec_id).get_best_peakgroup().get_normalized_retentiontime()
    retention_time[refrun_id][i] = refRT
    # Get the retention time of the precursor from reference run
    # refRT = best_feature_RT()
    # Iterate through all other runs and align them to the reference run
    for eXprun in eXps:
        ## Extract XICs from experiment run
        eXprun_id = eXprun.get_id()
        chrom_ids = run_chromIndex_map[eXprun_id][prec_id]
        if invalid_chromIndex in chrom_ids:
            print("Can't extract XICs for precursor {} from run {}".format(prec_id, refrun_id))
            continue
        XICs_eXp = extractXIC_group(chrom_file_accessor[eXprun_id], chrom_ids)
        XICs_eXp_sm = smoothXICs(XICs_eXp, sm)
        # Get time component
        t_ref = XICs_ref_sm[0][0]
        t_eXp = XICs_eXp_sm[0][0]
        ## Set up constraints for penalizing similarity matrix
        # Get RT tranfromation from refrun to eXprun.
        globalAligner, stdErr = get_pair_trafo(tr_data, spl_aligner, eXprun, refrun, multipeptides, max_rt_diff = 300)
        # Get constraining end-points using Global alignment
        B1p = globalAligner.predict([ t_ref[0] ])[0]
        B2p = globalAligner.predict([ t_ref[-1] ])[0]
        # Get width of the constraining window
        adaptiveRT = RSEdistFactor * stdErr
        samplingTime = (t_ref[-1] - t_ref[0]) / (len(t_ref) - 1)
        noBeef = np.ceil(adaptiveRT/samplingTime)
        ## Get intensity values to build similarity matrix
        intensityList_ref = np.array([xic[1] for xic in XICs_ref_sm], dtype=np.double)
        intensityList_eXp = np.array([xic[1] for xic in XICs_eXp_sm], dtype=np.double)
        intensityList_ref = np.ascontiguousarray(intensityList_ref)
        intensityList_eXp = np.ascontiguousarray(intensityList_eXp)
        chromAlignObj = AffineAlignObj(256, 256)
        alignChromatogramsCpp(chromAlignObj, intensityList_ref, intensityList_eXp,
                                    alignType = b"hybrid", tA = t_ref, tB = t_eXp,
                                    normalization = b"mean", simType = b"dotProductMasked",
                                    B1p = B1p, B2p = B2p, noBeef = noBeef,
                                    goFactor = 0.125, geFactor = 40,
                                    cosAngleThresh = 0.3, OverlapAlignment = True,
                                    dotProdThresh = 0.96, gapQuantile = 0.5,
                                    hardConstrain = False, samples4gradient = 100)
        AlignedIndices = pd.DataFrame(list(zip(chromAlignObj.indexA_aligned, chromAlignObj.indexB_aligned)),
	                                    columns =['indexAligned_ref', 'indexAligned_eXp'])
        AlignedIndices['indexAligned_ref'] -= 1
        AlignedIndices['indexAligned_eXp'] -= 1
        t_ref_aligned = get_aligned_time(t_ref, AlignedIndices['indexAligned_ref'])
        t_ref_aligned = pd.Series(t_ref_aligned).interpolate(method ='linear', limit_area = 'inside')
        t_eXp_aligned = get_aligned_time(t_eXp, AlignedIndices['indexAligned_eXp'])
        t_eXp_aligned = pd.Series(t_eXp_aligned).interpolate(method ='linear', limit_area = 'inside')
        t_ref_aligned = np.asarray(t_ref_aligned)
        index = np.nanargmin(np.abs(t_ref_aligned - refRT))
        eXpRT = t_eXp_aligned[index]
        retention_time[eXprun_id][i] = eXpRT

retention_time = pd.DataFrame.from_dict(retention_time)
retention_time = retention_time.dropna(how='all')
"""
B1p = 19.955
B2p = 49.998
noBeef = 20.0
t_ref = np.array([9.9, 13.3, 16.7, 20.1, 23.5, 26.9, 30.4, 33.8, 37.2, 40.6])
t_eXp = np.array([9.9, 13.3, 16.7, 20.1, 23.5, 26.9, 30.4, 33.8, 37.2, 40.6])
intensityList_ref = np.array([[1,2,3,2,1,0,0,0,0,0]], dtype = np.double)
intensityList_eXp = np.array([[0,0,0,0,1,2,3,2,1,0]], dtype = np.double)
chromAlignObj = AffineAlignObj(256, 256)

chromAlignObj.indexA_aligned # [0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
chromAlignObj.indexB_aligned # [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 10]
"""

fig = plt.figure()
XIC_Group = XICs_ref_sm
for k in range(len(XIC_Group)):
    x = XIC_Group[k][0]
    y = XIC_Group[k][1]
    plt.plot(x, y)

plt.ylabel('Intensity')
plt.show()
plt.close()
