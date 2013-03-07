set -e
###################################
# Variables
best_score=d_score
write_all_pg=0 # 0 means no, 1 means yes
mpr_num_xval=5
export OMP_NUM_THREADS=8
export memusage=4096
export min_upper_edge=1
export strict="-no-strict"
sleeptime=60

echo "Min upper edge" $min_upper_edge
echo "file"  $file_path
echo "base"  $file_basename
echo "iRT"  $irt_library
echo "ini"  $ini
echo "codedir"  $codedir
echo "OpenMS"  $openms_dir
echo "num xval"  $mpr_num_xval
echo "fdr_cutoff"  $fdr_cutoff
echo "chrom_extraction_window"  $extraction_window


# count how many files / swathes we have to process
mytmp=`ls ${file_path}/${file_basename}*.mzML.gz | wc` 
nr_files=`echo $mytmp | cut -f 1 -d ' '`

# combine them
bsub -I -n $OMP_NUM_THREADS -R "rusage[mem=4096]" ${openms_dir}FileMerger -in ${file_basename}*_.featureXML -out ${file_basename}_combined.featureXML -threads $OMP_NUM_THREADS
tar czf runlogs_mrmpeakpicker.tar.gz lsf.o*
rm  lsf*
rm ${file_basename}*_.featureXML

echo "Converting to short format"

bsub -I -n $OMP_NUM_THREADS -R "rusage[mem=6144]" OpenSwathFeatureXMLToTSV -tr $library -in ${file_basename}_combined.featureXML -out ${file_basename}_combined.short_format.csv -short_format -threads $OMP_NUM_THREADS

# run mprophet and feed the results back into the featureXML with a specified FDR

echo "Running mProphet"
bsub -I -n $OMP_NUM_THREADS -W8:00 -R "rusage[mem=6144]"  R --slave --args bin_dir=/cluster/apps/openms/openswath-testing/mapdiv/scripts/mProphet/ mquest=${file_basename}_combined.short_format.csv workflow=LABEL_FREE num_xval=$mpr_num_xval run_log=FALSE write_classifier=0 write_all_pg=1 help=0 project=${file_basename}_mprophet < /cluster/apps/openms/openswath-testing/mapdiv/scripts/mProphet/mProphet.R > ${file_basename}_mprophet.mProphet
echo "Done Running mProphet"
mv ${file_basename}_mprophet_all_peakgroups.xls ${file_basename}_all_peakgroups.xls
tar czvf runlogs_mprophet.tar.gz ${file_basename}_mprophet*
rm ${file_basename}_mprophet*
#rm ${file_basename}_combined.short_format.csv

bsub -I -n $OMP_NUM_THREADS -W8:00 -R "rusage[mem=6144]" OpenSwathRewriteToFeatureXML -FDR_cutoff $fdr_cutoff -csv ${file_basename}_all_peakgroups.xls \
  -featureXML ${file_basename}_combined.featureXML -out ${file_basename}_combined_fdr.featureXML -threads $OMP_NUM_THREADS
bzip2 ${file_basename}_all_peakgroups.xls
bsub -I -n $OMP_NUM_THREADS OpenSwathFeatureXMLToTSV -in ${file_basename}_combined_fdr.featureXML -tr $library \
  -out ${file_basename}_combined_fdr.long_format.tsv -best_scoring_peptide $best_score -threads $OMP_NUM_THREADS
bsub -I -n $OMP_NUM_THREADS OpenSwathFeatureXMLToTSV -in ${file_basename}_combined_fdr.featureXML -tr $library \
  -out ${file_basename}_combined_fdr.short_format.tsv -best_scoring_peptide $best_score -short_format -threads $OMP_NUM_THREADS

##### Step 5
# look at the result 
uptime
echo "Done"
outname=${file_basename}_combined_fdr.short_format.tsv
for fdr in 0.20 0.15 0.10 0.05 0.01 0.002 0.001; do
  python /cluster/apps/openms/openswath-testing/mapdiv/scripts/fdr_cutoff.py $outname ${outname}_fdr.csv m_score DECOY_ $fdr FALSE TRUE 
  grep -v DECOY ${outname}_fdr.csv > ${outname}_fdr_nodc.csv
  python /cluster/apps/openms/openswath-testing/mapdiv/scripts/count_pep_prot.py ${outname}_fdr_nodc.csv
done;

##### Step 6
# clean up, zip and combine the featureXMLs and the chromatograms
bzip2 ${file_basename}_combined_fdr.featureXML
bzip2 ${file_basename}_combined.featureXML
tar cjf ${file_basename}_chromatograms.tar.bz2 *._chrom.mzML
bsub -I -n $OMP_NUM_THREADS -W8:00 ${openms_dir}FileMerger -in ${file_basename}*._chrom.mzML -out ${file_basename}_combined.chrom.mzML -threads $OMP_NUM_THREADS
bzip2 ${file_basename}_combined.chrom.mzML
rm ${file_basename}*._chrom.mzML 


