####################################################################################################
# run in a folder where code points to mapdiv code. E.g. code/ExtractChromatogram is executable

# you have to set the following environmental variables using export

# make sure that ${file_path}/${file_basename}*.mzML.gz contains all your swath files

# Sample usage

set -e
###################################
# Variables
best_score=d_score
write_all_pg=0 # 0 means no, 1 means yes
mpr_num_xval=5
export OMP_NUM_THREADS=4
export memusage=4096
export min_upper_edge=1
export strict="-no-strict"
sleeptime=60

# echo "Min upper edge" $min_upper_edge
echo "file"  $file_path
echo "base"  $file_basename
echo "iRT"  $irt_library
echo "ini"  $ini
echo "codedir"  $codedir
echo "OpenMS"  $openms_dir
echo "num xval"  $mpr_num_xval
echo "fdr_cutoff"  $fdr_cutoff
echo "mprophet-mode"  $mprophet

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
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${codedir}OpenSwathChromatogramExtractor -is_swath -in $f -tr $irt_library \
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
${codedir}OpenSwathRTNormalizer -in ${file_basename}.rtnorm.chrom.mzML -tr $irt_library -out ${file_basename}.rtnorm.trafoXML
tar czf runlogs_rt_normalizer.tar.gz lsf.o*
rm lsf.o*
rt_norm=${file_basename}.rtnorm.trafoXML

##### Step 1
# Run ExtractChromatogram
for f in ${file_path}/${file_basename}*.mzML.gz; do
  outname=`basename $f .mzML.gz`._chrom.mzML
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${codedir}OpenSwathChromatogramExtractor -is_swath -in $f -tr $library \
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
  bsub -W $maxtime -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${codedir}OpenSwathAnalyzer -in $infile -swath_files $f \
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

### ## # remove -inf statements
### ## # TODO remove nan statements
### ## for f in ${file_basename}*_.featureXML; do
### ##   sed -i 's/-inf/0/' $f
### ## done

# combine them
if [ "$mprophet" == "combined" ]; then
	source mProphet_combined.sh
else
	source mProphet_separate.sh
fi

