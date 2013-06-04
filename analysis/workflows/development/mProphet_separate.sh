
set -e
###################################
# Variables
best_score=d_score
write_all_pg=0 # 0 means no, 1 means yes
mpr_num_xval=5
export OMP_NUM_THREADS=4
export memusage=4096
# export min_upper_edge=1
export strict="-no-strict"
sleeptime=60

# echo "Min upper edge" $min_upper_edge
echo "file"  $file_path
echo "base"  $file_basename
# echo "iRT"  $irt_library
#echo "ini"  $ini
echo "codedir"  $codedir
echo "OpenMS"  $openms_dir
echo "num xval"  $mpr_num_xval
echo "fdr_cutoff"  $fdr_cutoff
# echo "chrom_extraction_window"  $extraction_window


# count how many files / swathes we have to process
mytmp=`ls ${file_path}/${file_basename}*.mzML.gz | wc` 
nr_files=`echo $mytmp | cut -f 1 -d ' '`

echo "Start OpenSWATH job for ${nr_files} files"
uptime

# combine them
#${openms_dir}/FileMerger -in ${file_basename}*_.featureXML -out ${file_basename}_combined.featureXML
tar czf runlogs_mrmpeakpicker.tar.gz lsf.o*
rm  lsf*
#rm lsf* #${file_basename}*_.featureXML

echo "Converting to short format"
uptime

##### Step 3 Convert to short format
## counter=0
for f in ${file_basename}*_.featureXML; do
  outname=`basename $f _.featureXML`_.short_format.csv
  #echo "do ${outname}"
  outname=`echo $outname | sed 's/[%~]//g'`
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${codedir}OpenSwathFeatureXMLToTSV -tr $library \
    -in $f -out ${outname} -short_format -threads $OMP_NUM_THREADS 2>&1 > /dev/null
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files are converted to TSV"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the conversion to TSV "
else 
    echo "Not all files successfully completed the conversion to TSV"
    exit
fi
uptime

tar czf runlogs_featurexmltotsv.tar.gz lsf.o*
rm lsf* #

##### Step 4 Run mProphet
echo "Run mprophet"
uptime

# If_Combined_fast
# combined fast is an option where we combine the csv files (instead of the
# featureXMLs) and then run mProphet on it. This is most likely much faster
# than combining featureXMLs since it takes over 30 minutes to combine the
# featureXMLs and at least as long to rewrite to the large, combined xml.
if [ "$mprophet" == "combined_fast" ]; then

echo "do mprophet combined fast"
outname=${file_basename}_combined.short_format.csv
for f in *_.short_format.csv; do
    head -n 1 $f > $outname
    break
done

for f in *_.short_format.csv; do
    tail -n +2 $f >> $outname
done

rm *_.short_format.csv


# run mprophet 
echo "Running mProphet"
bsub -I -n $OMP_NUM_THREADS -W8:00 -R "rusage[mem=6144]"  R --slave --args bin_dir=mProphet/ mquest=${file_basename}_combined.short_format.csv workflow=LABEL_FREE num_xval=$mpr_num_xval run_log=FALSE write_classifier=0 write_all_pg=1 help=0 project=${file_basename}_mprophet < mProphet/mProphet.R > ${file_basename}_mprophet.mProphet
echo "Done Running mProphet"
uptime
# mv ${file_basename}_mprophet_peakgroups.xls ${file_basename}_all_peakgroups.xls
tar czvf runlogs_mprophet.tar.gz ${file_basename}_mprophet*
rm ${file_basename}_mprophet*

outname=${file_basename}_all_peakgroups.xls

# If_Combined_fast
else
echo "do mprophet separate "

### ## ## Lets try mprophet mutliple times.. might just work
### ## for j in {0..5}; do

# run mprophet and feed the results back into the featureXML with a specified FDR
for f in *_.short_format.csv;
do
    if [ -s $f ]
    then 
    bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" "R --slave --args bin_dir=mProphet/ \
      mquest=${f} workflow=LABEL_FREE num_xval=$mpr_num_xval \
      run_log=FALSE write_classifier=1 write_all_pg=${write_all_pg} help=0 project=${f}_mprophet < \
      mProphet/mProphet.R > ${f}_mprophet.mProphet" \
      2>&1 > /dev/null;
    else
        echo File $f is empty, continue;
        echo "done"
	let nr_files=nr_files-1;
    fi
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files are run through mProphet"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

### ### files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
### ### files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
### ### if [ "$files_lsf" == "$nr_files" ]; then
### ###     echo "All files successfully completed the mProphet run"
### ###     break
### ### else 
### ###     echo "Not all files successfully completed the $j mProphet run, lets try again"
### ###     tar czf runlogs_lsf_mprophet_$j.tar.gz lsf.o*
### ###     rm lsf* 
### ### fi
### ### 
### ### done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the mProphet run"
else 
    echo "Not all files successfully completed the mProphet run, abort!"
    exit
fi

uptime

outname=result_allpeakgroups_${file_basename}_mProphet_combined.csv
for f in *.csv_mprophet_peakgroups.xls; do
    head -n 1 $f > $outname
    break
done

for f in *.csv_mprophet_peakgroups.xls; do
    tail -n +2 $f >> $outname
done

tar czf runlogs_lsf_mprophet.tar.gz lsf.o*
tar czf runlogs_mprophet.tar.gz *csv_mprophet*
rm *csv_mprophet*
rm lsf* 

# If_Combined_fast
fi

##### Step 5
# look at the result 
uptime
echo "Done"
# sed -i 's/\tNA$/\t100/' $outname
#for fdr in 0.20 0.15 0.10 0.05 0.025 0.01 0.005 0.002 0.001; do
for fdr in 0.20 0.15 0.10 0.05 0.01 0.002 0.001; do
  python fdr_cutoff.py $outname ${outname}_fdr.csv m_score DECOY_ $fdr FALSE TRUE 
  grep -v DECOY ${outname}_fdr.csv > ${outname}_fdr_nodc.csv
  python count_pep_prot.py ${outname}_fdr_nodc.csv
done;

python fdr_cutoff.py $outname ${outname}_fdr.csv m_score DECOY_ $fdr_cutoff FALSE TRUE 
grep -v DECOY ${outname}_fdr.csv > ${outname}_fdr_${fdr_cutoff}_nodecoy.csv


##### Step 6
# Feed the results back into the featureXML with a specified FDR
echo "Rewriting to FeatureXML"
uptime
mprophet_outfile=$outname
for f in ${file_basename}*_.featureXML; do
  outname=`basename $f _.featureXML`_fdr.featureXML
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=6144]" OpenSwathRewriteToFeatureXML -FDR_cutoff $fdr_cutoff -csv $mprophet_outfile \
    -featureXML $f -out ${outname} -threads $OMP_NUM_THREADS
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files were re-written to featureXML"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the re-writing to featureXML "
else 
    echo "Not all files successfully completed the re-writing to featureXML "
    exit
fi
uptime

tar czf runlogs_rewrite_toFeatureXML.tar.gz lsf.o*
rm lsf* 

echo "Merge all featureXML with fdr"
uptime
bsub -I -n $OMP_NUM_THREADS -R "rusage[mem=4096]" ${openms_dir}FileMerger -in ${file_basename}*_fdr.featureXML -out ${file_basename}_combined__fdr.featureXML -threads $OMP_NUM_THREADS
rm ${file_basename}*_fdr.featureXML
mv ${file_basename}_combined__fdr.featureXML ${file_basename}_combined_fdr.featureXML

##### Step 7
# clean up, zip and combine the featureXMLs and the chromatograms
echo "Will start to clean up and zip files"
uptime
tar cjf ${file_basename}_features.tar.bz2 *_.featureXML
rm *_.featureXML
tar cjf ${file_basename}chromatograms.tar.bz2 *._chrom.mzML
rm ${file_basename}*._chrom.mzML 


#### add up the runtimes
## for f in runlogs_*; do tar -O -xzvf $f 2>/dev/null | grep 'CPU time   :'  | sed -n 's/CPU time.*:\(.*\)sec.*/\1  /p' | awk '{sum += $1} END {print sum}'; done
## add them all, dont forget the lsf master
## check memory tar -O -xzvf  runlogs_mrmpeakpicker.tar.gz 2>/dev/null | grep 'Max Memory :'  | sort

