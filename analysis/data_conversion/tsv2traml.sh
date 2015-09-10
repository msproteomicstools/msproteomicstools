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
# --------------------------------------------------------------------------
# $Maintainer: Lorenz Blum$
# $Authors: Lorenz Blum$
# --------------------------------------------------------------------------
set -e

[ "$#" -lt 2 ] && echo "Usage: $0 input.tsv output.traML
Converts input.tsv transition list into SWATH ready output.traML with decoys" && exit 1

#check inputfile
if [ $(grep -c "	" $1) -eq 0 ]
then
	echo "ERROR: No tabs found in infile, please make sure its tsv (not csv)!"
	exit 1
fi

if [ $(head -1 $1 | grep -c Tr_recalibrated) -ne 1 ] 
then
	echo "ERROR: Bad infile, header Tr_recalibrated not found" 
	exit 1
fi

#make TSV2TRAML compatible (ensure .csv suffix and fix line endings)
tempcsv=$1.temp.csv
cp $1 $tempcsv
dos2unix -q $tempcsv
mac2unix -q $tempcsv

#convert to traml
temptraml=$1.temp.traML
ConvertTSVToTraML -no_progress -in $tempcsv -out $temptraml

if [ $(cat $temptraml | wc -l) -le 7 ] 
then 
	echo "ERROR: Empty outputfile, no transitions!" 
	exit 1
fi

#add decoys
OpenSwathDecoyGenerator -in $temptraml -out $2 -threads 2 -method shuffle -append -exclude_similar -remove_unannotated
echo
echo "CREATED $2"

#cleanup
rm $tempcsv $temptraml