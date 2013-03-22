#!/bin/bash
set -e

input=
outdir=.
threads=1
windows=32
noms1map=

while getopts i:t:o:w:n opt
do
    case $opt in
    i)	input=$OPTARG;;
	o)	outdir=$OPTARG;;
    t)  threads=$OPTARG;;
	n)	noms1map="noms1map";;
    ?)  echo "Usage: $0 -i input.mzXML [-o outdir] [-t threads] [-w numWnd] [-n] 
Splits input.mzXML into its windows and converts them into mzML.gz using FileConverter and gzip. 
[-o outdir] location where split.mzML.gz are written. default=current dir
[-n] prevents writing of ms1map
[-w numWnd] number of swath windows. default=32
[-t numThreads] parallelizes the process using multiple processors. default=1"
        exit 1;;
    esac
done

#split
split_mzXML_intoSwath.py $input $windows $outdir $noms1map

#convert
for i in $outdir/split*.mzXML
do
	echo ${i%%.mzXML}
done | xargs -P $threads -I file sh -c '{ FileConverter -in "file.mzXML" -out "file.mzML"; gzip -fv "file.mzML"; rm -v "file.mzXML"; }' 
