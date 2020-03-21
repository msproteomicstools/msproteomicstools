#!/usr/bin/env python
# -*- coding: utf-8  -*-
import os, sys, csv, time

import argparse
import msproteomicstoolslib.math.Smoothing as smoothing
from msproteomicstoolslib.version import __version__ as version
from msproteomicstoolslib.math.chauvenet import chauvenet

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
from msproteomicstoolslib.format.TransformationCollection import initialize_transformation
from pyopenms import OnDiscMSExperiment
from analysis.chromatogram_utils.smooth_chromatograms import chromSmoother
from analysis.chromatogram_utils.chromatogramMapper import *
from analysis.alignment.reference_run_selection import referenceForPrecursor
from DIAlignPy import AffineAlignObj, alignChromatogramsCpp
from msproteomicstoolslib.format.TransformationCollection import TransformationCollection, LightTransformationData
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import write_out_matrix_file, get_pair_trafo
from msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm import AlignmentAlgorithm

def get_aligned_time(t, indices, skipValue = -1):
    t_aligned = [None]*len(indices)
    for i in range(len(indices)):
        index = indices[i]
        if index != skipValue:
            t_aligned[i] = t[index]
    t_aligned =  np.asarray(pd.Series(t_aligned).interpolate(method ='linear', limit_area = 'inside'))
    return t_aligned

def updateRetentionTime(run, peptide_group_label, prec_id, t_ref_aligned, t_eXp_aligned):
    prec = run.getPrecursor(peptide_group_label, prec_id)
    if prec is None:
        # TODO
        # Map boundaries, integrate area.
        # Create a precursor with the new peak-group and add it to the run. 
        return
    else:
        # Get their retention times.
        mutable = [list(pg) for pg in prec.peakgroups_]
        for k in range(len(mutable)):
            # Get corresponding retention time from reference run.
            index = np.nanargmin(np.abs(t_eXp_aligned - mutable[k][2]))
            # Update the retention time.
            mutable[k][2] = t_ref_aligned[index]
        # Insert back in precusor.
        prec.peakgroups_ = [ tuple(m) for m in mutable]

def RTofAlignedXICs(XICs_ref_sm, XICs_eXp_sm, RSEdistFactor, tr_data, spl_aligner, eXprun, refrun, multipeptides,
                            alignType = b"hybrid", normalization = b"mean", simType = b"dotProductMasked",
                            goFactor = 0.125, geFactor = 40,
                            cosAngleThresh = 0.3, OverlapAlignment = True,
                            dotProdThresh = 0.96, gapQuantile = 0.5,
                            hardConstrain = False, samples4gradient = 100):
    # Get time component
    t_ref = XICs_ref_sm[0][0]
    t_eXp = XICs_eXp_sm[0][0]
    ## Set up constraints for penalizing similarity matrix
    # Get RT tranfromation from refrun to eXprun.
    globalAligner, stdErr = get_pair_trafo(tr_data, spl_aligner, eXprun, refrun, multipeptides)
    # Get constraining end-points using Global alignment
    B1p = globalAligner.predict([ t_ref[0] ])[0]
    B2p = globalAligner.predict([ t_ref[-1] ])[0]
    # Get width of the constraining window
    adaptiveRT = RSEdistFactor * stdErr
    samplingTime = (t_ref[-1] - t_ref[0]) / (len(t_ref) - 1)
    noBeef = np.ceil(adaptiveRT/samplingTime)
    ## Get intensity values to build similarity matrix
    intensityList_ref = np.ascontiguousarray(np.array([xic[1] for xic in XICs_ref_sm], dtype=np.double))
    intensityList_eXp = np.ascontiguousarray(np.array([xic[1] for xic in XICs_eXp_sm], dtype=np.double))
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
    t_eXp_aligned = get_aligned_time(t_eXp, AlignedIndices['indexAligned_eXp'])
    return t_ref_aligned, t_eXp_aligned

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
MZs = mzml_accessors(runs, MStoFeature)
MZs.get_precursor_to_chromID(precursor_to_transitionID)

this_exp = Experiment()
this_exp.set_runs(runs)
start = time.time()
fdr_cutoff = 0.05
multipeptides = this_exp.get_all_multipeptides(fdr_cutoff, verbose=False, verbosity= 10)
print("Mapping the precursors took %0.2fs" % (time.time() - start) )

# Reference based alignment
best_run = this_exp.determine_best_run(alignment_fdr_threshold = 0.05)
reference_run = referenceForPrecursor(refType="precursor_specific", run= best_run, alignment_fdr_threshold = 0.05).get_reference_for_precursors(multipeptides)
# Pairwise global alignment
spl_aligner = SplineAligner(alignment_fdr_threshold = 0.05, smoother="lowess", experiment=this_exp)
tr_data = initialize_transformation()
# Initialize XIC smoothing function
chrom_smoother = chromSmoother(smoother = "sgolay", kernelLen = 11, polyOrd = 4)

# TODO Use multiprocessing tool for extracting chromatograms
RSEdistFactor = 4
# Calculate the aligned retention time for each precursor across all runs
prec_ids = list(precursor_to_transitionID.keys())

for i in range(len(prec_ids)):
    prec_id = prec_ids[i] #9719 9720
    refrun = reference_run.get(prec_id)
    if not refrun:
        print("The precursor {} doesn't have any associated reference run. Skipping!".format(prec_id))
        continue
    eXps = list(set(runs) - set([refrun]))
    # Extract XICs from reference run and smooth it.
    refrun_id = refrun.get_id()
    XICs_ref = MZs.extractXIC_group(refrun, prec_id)
    if not XICs_ref:
        continue
    XICs_ref_sm = chrom_smoother.smoothXICs(XICs_ref)
    # For each precursor, we need peptide_group_label and trgr_id
    peptide_group_label = precursor_to_transitionID[prec_id][0]
    # Iterate through all other runs and align them to the reference run
    for eXprun in eXps:
        ## Extract XICs from experiment run
        eXprun_id = eXprun.get_id()
        XICs_eXp = MZs.extractXIC_group(eXprun, prec_id)
        if not XICs_eXp:
            continue
        XICs_eXp_sm = chrom_smoother.smoothXICs(XICs_eXp)
        t_ref_aligned, t_eXp_aligned = RTofAlignedXICs(XICs_ref_sm, XICs_eXp_sm, RSEdistFactor,
                            tr_data, spl_aligner, eXprun, refrun, multipeptides,
                            alignType = b"hybrid", normalization = b"mean", simType = b"dotProductMasked",
                            goFactor = 0.125, geFactor = 40,
                            cosAngleThresh = 0.3, OverlapAlignment = True,
                            dotProdThresh = 0.96, gapQuantile = 0.5,
                            hardConstrain = False, samples4gradient = 100)
        # Update retention time of all peak-groups to reference peak-group
        updateRetentionTime(eXprun, peptide_group_label, prec_id, t_ref_aligned, t_eXp_aligned)

AlignmentAlgorithm().align_features(multipeptides, rt_diff_cutoff = 40, fdr_cutoff = 0.01,
                    aligned_fdr_cutoff = 0.05, method = "best_overall")
al = this_exp.print_stats(multipeptides, 0.05, 0.1, 1)

fraction_needed_selected = 0.0
write_out_matrix_file("GDYD.csv", runs, multipeptides, fraction_needed_selected, "RT", 0.05)

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


"""
retention_time
    precursor_id  125704171604355508  6752973645981403097  2234664662238281994
0           1967             5050.60               5049.6              5049.60
1           2474             6458.33               Nan                 6458.60
2           3864             3701.50               3664.0              3701.29
3           4618             5222.12               5220.8              5220.80
4           9719             2794.40               2793.6              2793.60
5           9720             2541.83               2540.2              2540.20
6          13678             4522.31               4521.3              4521.30
7          13679             3874.85               3876.1              3876.10
8          14378             4305.99               4306.5              4306.50
9          17186             4649.70               4651.0              4651.00
10         20003             4131.70               4128.3              4130.72

prec_id = 3864 # 2474, 20003
# For 3864 picks a peak from run1 based on FDR, if matched over boundary the peak should be different that is
# not even picked by OpenSWATH. Same for run0, OpenSWATH didn't even pick that feature.
# For 2474, some reason it didn't even have reasonable peak-groups in all three runs. Investigate
runs[0].getPrecursor(precursor_to_transitionID[prec_id][0], prec_id).peakgroups_
runs[1].getPrecursor(precursor_to_transitionID[prec_id][0], prec_id).peakgroups_
runs[2].getPrecursor(precursor_to_transitionID[prec_id][0], prec_id).peakgroups_

run_id = runs[2].get_id()
chrom_ids = run_chromIndex_map[run_id][prec_id]
XICs_eXp = extractXIC_group(chrom_file_accessor[run_id], chrom_ids)
XICs_eXp_sm = chrom_smoother.smoothXICs(XICs_eXp)
fig = plt.figure()
XIC_Group = XICs_eXp_sm
for k in range(len(XIC_Group)):
    x = XIC_Group[k][0]
    y = XIC_Group[k][1]
    plt.plot(x, y)

plt.ylabel('Intensity')
plt.show()
plt.close()
"""