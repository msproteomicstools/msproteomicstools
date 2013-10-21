#!/bin/bash
# -*- coding: utf-8  -*-
# =========================================================================
#         msproteomicstools -- Mass Spectrometry Proteomics Tools
# =========================================================================
#
# Copyright (c) 2013, ETH Zurich
# For a full list of authors, refer to the file AUTHORS.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#2CPU  28G 53m
#4CPU  24G 32m
#8CPU  29G 31m
#16CPU 43G 15m

# --------------------------------------------------------------------------
# $Maintainer: Lorenz Blum$
# $Authors: Lorenz Blum$
# --------------------------------------------------------------------------

set -e

input=
outdir=.
threads=1
windows=32
noms1map=

while getopts i:t:o:w:n opt
do
    case $opt in
    i)  input=$OPTARG;;
    o)  outdir=$OPTARG;;
    t)  threads=$OPTARG;;
    n)  noms1map="noms1map";;
    w)  windows=$OPTARG;;
    ?)  echo "Usage: $0 -i input.mzXML [-o outdir] [-t threads] [-w numSwathes] [-n]
Splits input.mzXML into its windows and converts them into mzML.gz using FileConverter and gzip.
[-o outdir] location where split.mzML.gz are written. default=current dir
[-n] prevents writing of ms1map
[-w numSwathes] number of swathes. default=32
[-t numThreads] parallelizes the process using multiple processors. default=1
Note: If \$TMPDIR is set it is used as temporary directory"
        exit 1;;
    esac
done


if [ -z "$TMPDIR" ]; then
    echo TMPDIR not set, using OUTDIR $outdir
    TMPDIR=$outdir
fi

#split
split_mzXML_intoSwath.py $input $windows $TMPDIR $noms1map

#convert in parallel using xargs -P
for i in $TMPDIR/split*.mzXML
do
    echo ${i%%.mzXML}
done | xargs -P $threads -I file sh -c '{ FileConverter -no_progress -in "file.mzXML" -out "file.mzML"; gzip -fv "file.mzML"; rm -v "file.mzXML"; }'

if [ "$TMPDIR" != "$outdir" ]; then
    echo Moving result files to outdir
    mv -v $TMPDIR/*.mzML.gz $outdir/
fi