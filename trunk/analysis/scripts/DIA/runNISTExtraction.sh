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
# $Maintainer: Hannes Roest $
# $Authors: Hannes Roest $
# --------------------------------------------------------------------------

# Shell script to convert parts of a NIST library format into TraML and use it to
# extract and score chromatograms using OpenSWATH.
# 
# 1. Downloaded NIST library data from http://peptide.nist.gov (manually)
# 2. set the MS2_location variable
# 3. run.sh input.txt extraction_window_size
# 

MS2_location="/media/data/tmp/MSE_2.mzML"

echo $1
echo $2
echo $MS2_location

python NISTParser.py $1.txt $1.csv
ConvertTSVToTraML -in $1.csv -out $1.TraML
OpenSwathChromatogramExtractor -tr $1.TraML -in $MS2_location -out $1_MS2.mzML -extraction_window $2
# OpenSwathChromatogramExtractor -tr $1.TraML -in /media/data/tmp/MSE_1.mzML -out $1_MS1.mzML -extraction_window $2

OpenSwathAnalyzer -in $1_corr_MS2.mzML -tr $1_corr.TraML -algorithm:Scores:use_rt_score false -out $1_corr.featureXML
OpenSwathFeatureXMLToTSV -in $1_corr.featureXML -tr $1_corr.TraML -out $1_corr_result.csv -short_format

