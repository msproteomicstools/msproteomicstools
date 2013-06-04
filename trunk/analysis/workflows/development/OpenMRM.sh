#!/usr/bin/env bash

# This assumes we are in a directory with n mzXML files and m TraML files that correspond to the mzXML files
# you have to set the following parameters before you can run it:
<< PARAMETERS
# If there are chromatograms that cannot be mapped (e.g. iRT chromatograms), one has to use no-strict
strict="-no-strict"
analyzer_ini="analyzer.ini"
mrmmapper_ini="mrmmapper.ini"
write_all_pg=0 # 0 means no, 1 means yes
mpr_num_xval=5
# for "by_name" to work, for each .mzXML file there needs to be an .TraML file with the exact same name (except the ending)
mapping="by_name" # by_name, "all" Mapping can be done by filename or all files can be combined
scoring="none" # can be none, ddb, separate, combined => no mProphet, run on each file or run on all files combined
fdr_cutoff=1.0
# if scoring is ddb, an assay library has to be provided
ddb_library=combined.TraML # if scoring is ddb, an assay library has to be provided
PARAMETERS

# module load imsbtools
# module load openms/svn
set -e

# Step 1: convert the mzXML into mzML files
for f in *.mzXML; 
  do
  bname=`basename $f .mzXML`
  FileConverter -in $f -out $bname.orig.mzML
done

# Step 2: convert the individual TraML files into one
FileMerger -in *.TraML -out combined.TraML

# Step 3: map the TraML to the mzML
for f in *.orig.mzML; 
  do
  bname=`basename $f .orig.mzML`
  if [ "$mapping" == "by_name" ]; then
    MRMMapper -in $f -out ${bname}_m.mzML -tr $bname.TraML -ini $mrmmapper_ini $strict 
  else
    MRMMapper -in $f -out ${bname}_m.mzML -tr combined.TraML -ini $mrmmapper_ini $strict
  fi
done

# Step 3: map for the iRT peptides 
for f in *.orig.mzML; 
  do
  bname=`basename $f .orig.mzML`
  MRMMapper -in $f -out ${bname}_mappedRT.mzML -tr $irt_library -ini $mrmmapper_ini $strict
done

# Step 4: do the RT normalisation
for f in *_mappedRT.mzML; 
  do
  bname=`basename $f _mappedRT.mzML`
  OpenSwathRTNormalizer -in $f -tr $irt_library -out $bname.trafoXML
done

# Step 5: do the peak picking
for f in *_m.mzML; 
  do
  bname=`basename $f _m.mzML`
  OpenSwathAnalyzer -in $f -tr combined.TraML -rt_norm $bname.trafoXML -out $bname.featureXML -ini $analyzer_ini -no-strict
done

if [ "$scoring" == "combined" ]; then

  # Step 6: merge and convert to TSV
  FileMerger -in *.featureXML -out combined.featureXML
  OpenSwathFeatureXMLToTSV -in combined.featureXML -out combined.csv -tr combined.TraML -short_format

  # Step 7: run mProphte and write back to featureXML 
  R --slave --args bin_dir=/cluster/apps/openms/openswath-master/OpenSWATH/scripts/mProphet/  \
      mquest=combined.csv workflow=LABEL_FREE num_xval=$mpr_num_xval \
      run_log=FALSE write_classifier=1 write_all_pg=${write_all_pg} help=0 project=combined_mprophet < \
      /cluster/apps/openms/openswath-master/OpenSWATH/scripts/mProphet/mProphet.R > log_mprophet.mProphet

  # TODO test!
  rewriteToFeatureXML -FDR_cutoff $fdr_cutoff -csv combined_mprophet_all_peakgroups.xls \
    -featureXML combined.featureXML -out combined_fdr.featureXML 

else
if [ "$scoring" == "separate" ]; then
    # separate mprophet

    # Step 6: convert to TSV
    for f in *.featureXML; 
      do
      bname=`basename $f .featureXML`
      OpenSwathFeatureXMLToTSV -in $f -out $bname.csv -tr combined.TraML -short_format
    done

  # Step 7: run mProphte and write back to featureXML 
    for f in *.csv;
    do
        R --slave --args bin_dir=/cluster/scratch/malars/openswath/mQuest_mProphet/libs/r_libs \
          mquest=${f} workflow=LABEL_FREE num_xval=$mpr_num_xval \
          run_log=FALSE write_classifier=1 write_all_pg=${write_all_pg} help=0 project=${f}_mprophet < \
          /cluster/scratch/malars/openswath/mQuest_mProphet/mQuest/bin/mProphet.R > ${f}_mprophet.mProphet
    done
    for f in *.csv;
    do
      bname=`basename $f .csv`
      rewriteToFeatureXML -FDR_cutoff $fdr_cutoff -csv ${file_basename}_all_peakgroups.xls \
        -featureXML ${file_basename}_combined.featureXML -out ${file_basename}_combined_fdr.featureXML 
    #rm ${file_basename}_all_peakgroups.xls
    done;

    # Step 8: create new, final CSV
    FeatureXMLToTSV -in ${file_basename}_combined_fdr.featureXML -tr $library \
      -out ${file_basename}_combined_fdr.long_format.tsv -best_scoring_peptide $best_score

else
if [ "$scoring" == "ddb" ]; then

    for f in *.featureXML; 
      do 
      echo $f;
      bname=`basename $f .featureXML`;
      OpenSwathConfidenceScoring  \
        -in  $f \
        -lib $ddb_library \
        -out  $bname.score.featureXML \
        -trafo $bname.trafoXML \
        # -ini 
      done

    for f in *.score.featureXML; 
    do
      bname=`basename $f score.featureXML`;
      OpenSwathFeatureXMLToTSV -in $f -out $bname.short_format.csv -tr combined.TraML -short_format
    done

else
    # no mProphet or ddb, just give out the textfiles
    # Step 6: convert to TSV
    for f in *.featureXML; 
      do
      bname=`basename $f .featureXML`
      OpenSwathFeatureXMLToTSV -in $f -out $bname.short_format.csv -tr combined.TraML -short_format
      OpenSwathFeatureXMLToTSV -in $f -out $bname.long_format.csv -tr combined.TraML 
    done
fi
fi
fi


