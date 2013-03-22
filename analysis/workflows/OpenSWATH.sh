# #########################################
# ########## ANALYSIS PARAMETERS ##########
# # Path to SWATH MS files:
# export file_path=/cluster/scratch_xl/shareholder/imsb_ra/openswath/data/mTB_run/olgas_L120527_011_SW/
# # SWATH MS files basename:
# export file_basename=split_olgas_L120527_011
# # Peptide assay library:
# export library=/cluster/scratch_xl/shareholder/imsb_ra/openswath/tramlpile/SamBader_Mtb_110_shuffle.TraML
# # iRT assay library (Do not change if iRT-kit was used):
# export irt_library=/cluster/scratch_xl/shareholder/imsb_ra/openswath/tramlpile/hroest_DIA_iRT.TraML
#  # FDR Cutoff:
# export fdr_cutoff=0.01

# #########################################
# ### EXPERT PARAMETERS - DO NOT CHANGE ###
# module load openms/openswath-loblum
# # Path to OpenSWATH.ini file:
# export ini=/cluster/apps/guse/20130215/msproteomictools/trunk/analysis/workflows/OpenSWATH.ini
# # Overlap of swaths in Da:
# export min_upper_edge=1
# # MS2 extraction window in Da:
# export extraction_window=0.05
# # RT extraction window in ±s (set to -1 to extract whole chromatogram):
# export rt_extraction_window=300
# # mProphet: specify main score:
# export main_var=xx_swath_prelim_score
# # mProphet: specify additional scores:
# export vars="bseries_score elution_model_fit_score intensity_score isotope_correlation_score isotope_overlap_score library_corr library_rmsd log_sn_score massdev_score massdev_score_weighted norm_rt_score xcorr_coelution xcorr_coelution_weighted xcorr_shape xcorr_shape_weighted yseries_score"
# # mProphet: specify LDA model (leave empty for none):
# export mprophet_model=""
# # mProphet: write all peak groups (0 means no, 1 means yes):
# export write_all_pg=0
# # mProphet: number of cross validation runs:
# export mpr_num_xval=5
# #########################################

# bsub -W 8:00 -R "rusage[mem=2048]" OpenSWATH.sh

set -e
###################################
# Variables

sleeptime=60

export strict="-no-strict"
export OMP_NUM_THREADS=4
export memusage=4096
export best_score=d_score

echo "Base configuration:"
echo "-------------------"
echo "OpenMS dir"  $openms_dir
echo "mProphet dir"  $mprophet_dir
echo "Number of threads"  $OMP_NUM_THREADS
echo "Amount of RAM per core"  $memusage
echo ""

echo "Run configuration:"
echo "------------------"
echo "Path to SWATH MS files:"  $file_path
echo "SWATH MS files basename:"  $file_basename
echo "Peptide assay library:"  $library
echo "iRT assay library:"  $irt_library
echo "FDR Cutoff:"  $fdr_cutoff
echo ""

echo "Analysis configuration:"
echo "-----------------------"
echo "Path to OpenSWATH.ini file:"  $ini
echo "Overlap of swaths in Da:" $min_upper_edge
echo "MS2 extraction window in Da:" $extraction_window
echo "RT extraction window in ±s:" $rt_extraction_window
echo "mProphet: specify main score:"  $main_var
echo "mProphet: specify additional scores:"  $vars
echo "mProphet: specify LDA model:"  $mprophet_model
echo "mProphet: write all peak groups:"  $write_all_pg
echo "mProphet: number of cross validation runs:"  $mpr_num_xval
echo ""

if [ -z "$maxtime" ];
then
    maxtime=59
fi;

# count how many files / swathes we have to process
mytmp=`ls ${file_path}/${file_basename}*.mzML.gz | wc` 
nr_files=`echo $mytmp | cut -f 1 -d ' '`

echo "Start OpenSWATH job for ${nr_files} files"
uptime

##### Step 0
# Run ExtractChromatogram for the RT Peptides
for f in ${file_path}/${file_basename}*.mzML.gz; do
  outname=`basename $f .mzML.gz`._rtnorm.chrom.mzML
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${openms_dir}OpenSwathChromatogramExtractor -is_swath -in $f -tr $irt_library \
    -out $outname -ini $ini -threads $OMP_NUM_THREADS -rt_extraction_window -1 2>&1 > /dev/null;
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files are iRT extracted"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the iRT extraction"
else
    echo "Not all files successfully completed the iRT extraction"
    exit
fi
uptime

# merge files and run RT Normalizer
${openms_dir}FileMerger -in ${file_basename}*._rtnorm.chrom.mzML -out ${file_basename}.rtnorm.chrom.mzML
rm ${file_basename}*._rtnorm.chrom.mzML 
${openms_dir}OpenSwathRTNormalizer -in ${file_basename}.rtnorm.chrom.mzML -tr $irt_library -out ${file_basename}.rtnorm.trafoXML
tar czf runlogs_rt_normalizer.tar.gz lsf.o*
rm lsf.o*
rt_norm=${file_basename}.rtnorm.trafoXML

##### Step 1
# Run ExtractChromatogram
for f in ${file_path}/${file_basename}*.mzML.gz; do
  outname=`basename $f .mzML.gz`._chrom.mzML
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${openms_dir}OpenSwathChromatogramExtractor -is_swath -in $f -tr $library \
    -out $outname -ini $ini -rt_norm $rt_norm -threads $OMP_NUM_THREADS 2>&1 > /dev/null;
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files are extracted"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the extraction"
else 
    echo "Not all files successfully completed the extraction"
    exit
fi
uptime

#${openms_dir}/FileMerger -in ${file_basename}*._chrom.mzML -out ${file_basename}_combined.chrom.mzML
tar czf runlogs_extract_chromatogram.tar.gz lsf.o*
#rm ${file_basename}*._chrom.mzML lsf.o*
rm lsf.o*

##### Step 2
for f in ${file_path}/${file_basename}*.mzML.gz; do
  outname=`basename $f .mzML.gz`_.featureXML
  infile=`basename $f .mzML.gz`._chrom.mzML
  bsub -W $maxtime -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${openms_dir}OpenSwathAnalyzer -in $infile -swath_files $f \
    -tr $library -out $outname -ini $ini -rt_norm $rt_norm -threads $OMP_NUM_THREADS $strict 2>&1 > /dev/null
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files are picked"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the feature finding"
else 
    echo "Not all files successfully completed the feature finding"
    exit
fi

uptime

# count how many files / swathes we have to process
mytmp=`ls ${file_path}/${file_basename}*.mzML.gz | wc` 
nr_files=`echo $mytmp | cut -f 1 -d ' '`

# combine them
bsub -I -n $OMP_NUM_THREADS -R "rusage[mem=4096]" ${openms_dir}FileMerger -in ${file_basename}*_.featureXML -out ${file_basename}_combined.featureXML -threads $OMP_NUM_THREADS
tar czf runlogs_mrmpeakpicker.tar.gz lsf.o*
rm  lsf*
rm ${file_basename}*_.featureXML

##### Step 3
# Conversion

echo "Converting to short format"

bsub -I -n $OMP_NUM_THREADS -R "rusage[mem=6144]" OpenSwathFeatureXMLToTSV -tr $library -in ${file_basename}_combined.featureXML -out ${file_basename}_combined.short_format.csv -short_format -threads $OMP_NUM_THREADS

##### Step 4
#mProphet

# select the scores for mProphet
bsub -I -n $OMP_NUM_THREADS -R "rusage[mem=4096]" ${tools_dir}mProphetScoreSelector.sh ${file_basename}_combined.short_format.csv ${main_var} ${vars}

# run mprophet and feed the results back into the featureXML with a specified FDR

echo "Running mProphet"
if [ "$mprophet_model" != ""]; then
  bsub -I -n $OMP_NUM_THREADS -W8:00 -R "rusage[mem=6144]" ${mprophet_dir}mProphetRunner.sh mquest=${file_basename}_combined.short_format.csv workflow=LABEL_FREE num_xval=$mpr_num_xval run_log=FALSE write_classifier=1 write_all_pg=${write_all_pg} help=0 project=${file_basename}_mprophet use_classifier=${mprophet_model} > ${file_basename}_mprophet.mProphet
else
  bsub -I -n $OMP_NUM_THREADS -W8:00 -R "rusage[mem=6144]" ${mprophet_dir}mProphetRunner.sh mquest=${file_basename}_combined.short_format.csv workflow=LABEL_FREE num_xval=$mpr_num_xval run_log=FALSE write_classifier=1 write_all_pg=${write_all_pg} help=0 project=${file_basename}_mprophet > ${file_basename}_mprophet.mProphet
fi

echo "Done Running mProphet"
mv ${file_basename}_mprophet_*peakgroups.xls ${file_basename}_all_peakgroups.xls
tar czvf runlogs_mprophet.tar.gz ${file_basename}_mprophet*
rm ${file_basename}_mprophet*
#rm ${file_basename}_combined.short_format.csv

bsub -I -n $OMP_NUM_THREADS -W8:00 -R "rusage[mem=6144]" OpenSwathRewriteToFeatureXML -FDR_cutoff $fdr_cutoff -csv ${file_basename}_all_peakgroups.xls \
  -featureXML ${file_basename}_combined.featureXML -out ${file_basename}_combined_fdr.featureXML -threads $OMP_NUM_THREADS
bzip2 ${file_basename}_all_peakgroups.xls
bsub -I -n $OMP_NUM_THREADS OpenSwathFeatureXMLToTSV -in ${file_basename}_combined_fdr.featureXML -tr $library \
  -out ${file_basename}_combined_fdr.long_format.csv -best_scoring_peptide $best_score -threads $OMP_NUM_THREADS
bsub -I -n $OMP_NUM_THREADS OpenSwathFeatureXMLToTSV -in ${file_basename}_combined_fdr.featureXML -tr $library \
  -out ${file_basename}_combined_fdr.short_format.csv -best_scoring_peptide $best_score -short_format -threads $OMP_NUM_THREADS

##### Step 5
# look at the result 
uptime
echo "Done"
outname=${file_basename}_combined_fdr.short_format.csv
for fdr in 0.20 0.15 0.10 0.05 0.01 0.002 0.001; do
 ${mprophet_dir}fdr_cutoff.py $outname ${outname}_fdr.csv m_score DECOY_ $fdr FALSE TRUE 
  grep -v DECOY ${outname}_fdr.csv > ${outname}_fdr_nodc.csv
  ${mprophet_dir}count_pep_prot.py ${outname}_fdr_nodc.csv
done;

##### Step 6
# clean up, zip and combine the featureXMLs and the chromatograms
bzip2 ${file_basename}_combined_fdr.featureXML
bzip2 ${file_basename}_combined.featureXML
tar cjf ${file_basename}_chromatograms.tar.bz2 *._chrom.mzML
bsub -I -n $OMP_NUM_THREADS -W8:00 ${openms_dir}FileMerger -in ${file_basename}*._chrom.mzML -out ${file_basename}_combined.chrom.mzML -threads $OMP_NUM_THREADS
bzip2 ${file_basename}_combined.chrom.mzML
rm ${file_basename}*._chrom.mzML 
