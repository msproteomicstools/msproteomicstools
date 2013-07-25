# This script can be used and adapted to truncate and convert the raw data.
cd /IMSB/users/georger/html/osw_peptides/data

cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_human/split_napedro_L120417_*/*_3.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_water/split_napedro_L120224_*/*_3.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_yeast/split_napedro_L120302_*/*_3.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_human/split_napedro_L120417_*/*_6.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_water/split_napedro_L120224_*/*_6.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_yeast/split_napedro_L120302_*/*_6.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_human/split_napedro_L120417_*/*_10.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_water/split_napedro_L120224_*/*_10.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_yeast/split_napedro_L120302_*/*_10.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_human/split_napedro_L120417_*/*_7.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_water/split_napedro_L120224_*/*_7.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_yeast/split_napedro_L120302_*/*_7.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_human/split_napedro_L120417_*/*_13.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_water/split_napedro_L120224_*/*_13.mzML.gz ./
cp /cluster/scratch_xl/shareholder/imsb_ra/openswath/data/AQUA_fixed_yeast/split_napedro_L120302_*/*_13.mzML.gz ./

for filename in *_3.mzML.gz
# 2273.999
do
   bsub FileFilter -in $filename -out "${filename%.*.*}_filtered.mzML.gz" -rt 2073:2473
done

for filename in *_6.mzML.gz
# 3821.858
do
   bsub FileFilter -in $filename -out "${filename%.*.*}_filtered.mzML.gz" -rt 3621:4021
done

for filename in *_10.mzML.gz
# 3530.576
do
   bsub FileFilter -in $filename -out "${filename%.*.*}_filtered.mzML.gz" -rt 3330:3730
done

for filename in *_7.mzML.gz
# 2517.863
do
   bsub FileFilter -in $filename -out "${filename%.*.*}_filtered.mzML.gz" -rt 2317:2717
done

for filename in *_13.mzML.gz
# 2321.417
do
   bsub FileFilter -in $filename -out "${filename%.*.*}_filtered.mzML.gz" -rt 2121:2521
done

for filename in *_filtered.mzML.gz
do
   bsub FileConverter -in $filename -out "${filename%.*.*}.dta2d"
done
