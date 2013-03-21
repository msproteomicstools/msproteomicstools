#!/bin/bash
set -e

remove=0
threads=1
while getopts rt: opt
do
    case $opt in
    r)    remove=1;;
    t)    threads=$OPTARG;;
    ?)    printf "Usage: $0 [-r] [-t int] mzXMLfiles...
Converts given mzXMLfiles to mzML.gz using FileConverter and gzip. 
-r if set original mzXML are removed (to save space)
-t int parallelizes process using int processors"
          exit 1;;
    esac
done
shift $(($OPTIND - 1))

echo Using $threads threads
for i in $@
do
	echo ${i%%.mzXML}
done | xargs -P $threads -I file sh -c '{ FileConverter -in "file.mzXML" -out "file.mzML"; gzip -fv "file.mzML"; }' 

[ $remove -eq 1 ] && rm -v $@
