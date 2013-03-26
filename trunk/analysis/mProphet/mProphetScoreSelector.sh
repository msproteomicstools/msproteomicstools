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

[ "$#" -eq 0 ] && echo "Usage: $0 file main_var [vars...]" && exit

file=$1
#remove previous var_ annotations, add main_var to sedstring
sedstring="1s/main_var_//;1s/var_//g;1s/$2/main_var_$2/;"
shift 2
for i in $@
do
    #add vars to sedstring
	sedstring="${sedstring}1s/${i}/var_${i}/;"
done
#execute sed command once only
sed -i"" "$sedstring" $file
