module load openms/svn-current

# This script can be used and adapted to truncate and convert the raw data.

cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120417_*/*_3._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120224_*/*_3._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120302_*/*_3._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120417_*/*_6._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120224_*/*_6._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120302_*/*_6._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120417_*/*_10._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120224_*/*_10._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120302_*/*_10._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120417_*/*_7._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120224_*/*_7._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120302_*/*_7._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120417_*/*_13._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120224_*/*_13._chrom.mzML ./
cp /cluster/scratch_xl/shareholder/biol/georger/osw/test/AQUASky_ShotgunLibrary_3t_345_sh/split_napedro_L120302_*/*_13._chrom.mzML ./# This script can be used and adapted to truncate and convert the raw data.

mkdir VGDTVLYGK
for f in *_3._chrom.mzML; do
  bsub "python ../filterChrom.py $f VGDTVLYGK/$f '^[^D].*_VGDTVLYGK.UniMod'; FileConverter -in ./VGDTVLYGK/$f -out VGDTVLYGK/${f}.dta2d"
done;

mkdir IADIQLEGLR
for f in *_6._chrom.mzML; do
  bsub "python ../filterChrom.py $f IADIQLEGLR/$f '^[^D].*IADIQLEGLR.UniMod'; FileConverter -in ./IADIQLEGLR/$f -out IADIQLEGLR/${f}.dta2d"
done;

mkdir TGGDEFDEAIIK
for f in *_10._chrom.mzML; do
  bsub "python ../filterChrom.py $f TGGDEFDEAIIK/$f '^[^D].*TGGDEFDEAIIK.UniMod'; FileConverter -in ./TGGDEFDEAIIK/$f -out TGGDEFDEAIIK/${f}.dta2d"
done;

mkdir LITVEGPDGAGK
for f in *_7._chrom.mzML; do
  bsub "python ../filterChrom.py $f LITVEGPDGAGK/$f '^[^D].*LITVEGPDGAGK.UniMod'; FileConverter -in ./LITVEGPDGAGK/$f -out LITVEGPDGAGK/${f}.dta2d"
done;

mkdir LVDEEGNDVTPEK
for f in *_13._chrom.mzML; do
  bsub "python ../filterChrom.py $f LVDEEGNDVTPEK/$f '^[^D].*LVDEEGNDVTPEK.UniMod'; FileConverter -in ./LVDEEGNDVTPEK/$f -out LVDEEGNDVTPEK/${f}.dta2d"
done;
