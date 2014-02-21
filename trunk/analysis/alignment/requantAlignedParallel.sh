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

requantArgs=""
for i in $@
do
    if [ $i == "-h" ]; then
        echo "About: Wrapper around requantAlignedValues to execute it in parallel.

Usage: Call this wrapper like an all-files-at-once requantAlignedValues run. It will automatically do parallelisation in the background.
If [--threads n] argument is added, this value instead of #CPU threads are used."
        exit 1
    fi
    #parse args to variables: --in => $in
	if [[ $i == --* ]]; then
		varToUse=${i#--}
	else
		declare "$varToUse=${!varToUse}$i "
	fi
    #don't forward out/do_single_run args because they're overwritten by this script
    if [ "$varToUse" != "out" -a "$varToUse" != "out_matrix" -a "$varToUse" != "do_single_run" -a "$varToUse" != "threads" ]; then
        requantArgs+="$i "
    fi
done

[ -z "$threads" ] && threads=$(nproc --all|| echo 1)
echo Using $threads threads

uniqtmpdir=$(mktemp -d)
for i in $in
do
	tempout=$(mktemp --tmpdir=$uniqtmpdir)
	gz="$(find $(dirname $i) -name \*.gz)"
	if [ "$gz" != "" ]; then
	    echo -n "gunzip -c $gz > ${gz%%.gz} && "
	fi
	echo -n requantAlignedValues.py $requantArgs --out $tempout --do_single_run $i
	if [ ! -z "$gz" ]; then
	    echo -n " && rm ${gz%%.gz}"
	fi
	echo
done | parallel --halt 2 -j $threads

awk "NR==1 || FNR!=1" $uniqtmpdir/* $peakgroups_infile > $out
rm -r $uniqtmpdir

if [ ! -z "$out_matrix" ]; then
    compute_full_matrix.py --in $out --out_matrix $out_matrix
fi