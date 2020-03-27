#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta reference_run_selection.py>
Copyright (C) 2020 Shubham Gupta

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

--------------------------------------------------------------------------
$Maintainer: Shubham Gupta$
$Authors: Shubham Gupta$
--------------------------------------------------------------------------
"""
import os, sys, csv, time
import numpy as np
import pandas as pd
import argparse
from msproteomicstoolslib.version import __version__ as version
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

def RTofAlignedXICs(XICs_ref_sm, XICs_eXp_sm, tr_data, spl_aligner, eXprun, refrun, multipeptides,
                            RSEdistFactor = 4, alignType = b"hybrid", normalization = b"mean", simType = b"dotProductMasked",
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

def plotXICs(XIC_Group):
    fig = plt.figure()
    for k in range(len(XIC_Group)):
        x = XIC_Group[k][0]
        y = XIC_Group[k][1]
        plt.plot(x, y)
    plt.ylabel('Intensity')
    plt.show()

def handle_args():
    usage = ""
    usage += "\nFeature Alignment -- version %d.%d.%d" % version
    usage += "\n"
    usage += "\nThis program aligns MS2 chromatograms to map peakgroups across runs. Post-alignment peakgroups below a certain FDR cutoff will be selcted."
    usage += "\n"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--chromatograms', dest="chromatogram_files", required=True, nargs = '+', help = 'A list of chromatogram files', metavar="file.chrom.mzML")
    parser.add_argument('--features', dest="feature_files", required=True, nargs = '+', help = 'A list of pyProphet output files containing all peakgroups (use quotes around the filenames)', metavar="file.osw")
    parser.add_argument('--file_format', default='osw', help="Input file format (osw)", metavar="")
    parser.add_argument("--out_matrix", dest="matrix_outfile", default="", help="Matrix containing one peak group per row (supports .csv, .tsv or .xlsx)", metavar="")
    parser.add_argument("--fdr_cutoff", dest="fdr_cutoff", default=0.01, type=float, help="Fixed FDR cutoff used for seeding (only assays where at least one peakgroup in one run is below this cutoff will be included in the result), see also target_fdr for a non-fixed cutoff", metavar='0.01')
    parser.add_argument("--min_fdr_quality", dest="nonaligned_fdr_cutoff", default=0.01, type=float, help="m-score cutoff below that peakgroups need not be evaluated with retention time alignment.", metavar='0.01')
    parser.add_argument("--max_fdr_quality", dest="aligned_fdr_cutoff", default=0.05, help="Extension m-score score cutoff, peakgroups of this quality will still be considered for alignment during extension", metavar='0.05')
    parser.add_argument("--max_rt_diff", dest="rt_diff_cutoff", default=30, help="Maximal difference in RT (in seconds) for two aligned features", metavar='30')
    parser.add_argument("--frac_selected", dest="min_frac_selected", default=0.0, type=float, help="Do not write peakgroup if selected in less than this fraction of runs (range 0 to 1)", metavar='0')
    parser.add_argument('--global_align_method', default='loess', help="Method to use for a global RT alignment of runs (linear, loess).")
    parser.add_argument('--chromatogram_align_method', default='hybrid', help="Method to use for RT alignment of chromatograms (local, global, hybrid).")
    parser.add_argument('--method', default='best_overall', help="Method to use for the clustering (best_overall).")
    parser.add_argument("--verbosity", default=0, type=int, help="Verbosity (0 = little)", metavar='0')
    parser.add_argument("--matrix_output_method", dest="matrix_output_method", default='RT', help="Which columns are written besides Intensity (none, RT, score, source or full)", metavar="")
    parser.add_argument('--force', action='store_true', default=False, help="Force alignment")
    parser.add_argument("--version", dest="version", default="", help="Print version and exit")

    if any([a.startswith("--version") for a in sys.argv]):
        print(usage)
        sys.exit()
    
    args = parser.parse_args(sys.argv[1:])

    if args.matrix_outfile == "":
        args.matrix_outfile = "DIAlignPy.csv"
    if args.min_frac_selected < 0.0 or args.min_frac_selected > 1.0:
        raise Exception("Argument frac_selected needs to be a number between 0 and 1.0")

    return args

def main(options):
    infiles = options.feature_files
    chromatograms = options.chromatogram_files

    readfilter = ReadFilter()
    file_format = 'openswath'
    readmethod = "minimal"

    reader = SWATHScoringReader.newReader(infiles, file_format, readmethod, readfilter,
                                            enable_isotopic_grouping = False, read_cluster_id = False)
    reader.map_infiles_chromfiles(chromatograms)
    runs = reader.parse_files()
    MStoFeature = MSfileRunMapping(chromatograms, runs)
    precursor_to_transitionID, precursor_sequence = getPrecursorTransitionMapping(infiles[0])
    MZs = mzml_accessors(runs, MStoFeature)
    MZs.set_precursor_to_chromID(precursor_to_transitionID)

    this_exp = Experiment()
    this_exp.set_runs(runs)
    start = time.time()
    fdr_cutoff = options.aligned_fdr_cutoff
    multipeptides = this_exp.get_all_multipeptides(fdr_cutoff, verbose=False, verbosity= 10)
    print("Mapping the precursors took %0.2fs" % (time.time() - start) )

    # Reference based alignment
    # best_run = this_exp.determine_best_run(alignment_fdr_threshold = 0.05)
    reference_run = referenceForPrecursor(refType="precursor_specific", alignment_fdr_threshold = options.fdr_cutoff).get_reference_for_precursors(multipeptides)
    # Pairwise global alignment
    spl_aligner = SplineAligner(alignment_fdr_threshold = fdr_cutoff, smoother="lowess", experiment=this_exp)
    tr_data = initialize_transformation()
    # Initialize XIC smoothing function
    chrom_smoother = chromSmoother(smoother = "sgolay", kernelLen = 11, polyOrd = 4)

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
        XICs_ref = MZs.extractXIC_group(refrun, prec_id)
        if not XICs_ref:
            continue
        XICs_ref_sm = chrom_smoother.smoothXICs(XICs_ref)
        # For each precursor, we need peptide_group_label and trgr_id
        peptide_group_label = precursor_sequence[prec_id][0]
        # Iterate through all other runs and align them to the reference run
        for eXprun in eXps:
            ## Extract XICs from experiment run and smooth it.
            XICs_eXp = MZs.extractXIC_group(eXprun, prec_id)
            if not XICs_eXp:
                continue
            XICs_eXp_sm = chrom_smoother.smoothXICs(XICs_eXp)
            t_ref_aligned, t_eXp_aligned = RTofAlignedXICs(XICs_ref_sm, XICs_eXp_sm,
                                tr_data, spl_aligner, eXprun, refrun, multipeptides, RSEdistFactor = 4,
                                alignType = b"hybrid", normalization = b"mean", simType = b"dotProductMasked",
                                goFactor = 0.125, geFactor = 40,
                                cosAngleThresh = 0.3, OverlapAlignment = True,
                                dotProdThresh = 0.96, gapQuantile = 0.5,
                                hardConstrain = False, samples4gradient = 100)
            # Update retention time of all peak-groups to reference peak-group
            updateRetentionTime(eXprun, peptide_group_label, prec_id, t_ref_aligned, t_eXp_aligned)

    AlignmentAlgorithm().align_features(multipeptides, rt_diff_cutoff = 40, fdr_cutoff = 0.01,
                        aligned_fdr_cutoff = options.aligned_fdr_cutoff, method = options.method)
    al = this_exp.print_stats(multipeptides, 0.05, 0.1, 1)
    write_out_matrix_file(options.matrix_outfile, runs, multipeptides, options.min_frac_selected,
                         options.matrix_output_method, True, 0.05, precursor_sequence)

if __name__=="__main__":
    options = handle_args()
    main(options)

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
"""