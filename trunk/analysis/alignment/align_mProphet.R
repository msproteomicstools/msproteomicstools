#!/usr/bin/env Rscript
# -*- coding: utf-8  -*-
# 
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
# $Maintainer: Hannes Roest$
# $Authors: Hannes Roest$
# --------------------------------------------------------------------------
# 

# Script to align two runs based on id
# Hannes Roest, January 2013
# 
# The input is a csv file with the following column names:
#  - transition_group_id, peak_group_rank, RT, decoy, m_score (fdr score)
# ./align_mProphet.R infile1 infile2 0.0001

# possible use:
# export reffile=strep_align/Strep0_Repl1_R02/split_hroest_K120808_all_peakgroups.xls
# for f in strep_align/*/*all_peakgroups.xls; do echo $f; ./align_mProphet.R $reffile $f 0.0001; done

# This requires the following headers: transition_group_id, peak_group_rank, RT, decoy
unique_merge_id = "transition_group_id"

# we make a model where we predict y (r2) from x (r1). Thus run1 is our "master"
get_smoothing_spline <- function(r1, r2, m_score_cutoff, filename)
{
  best_r1 = subset(r1, m_score < m_score_cutoff)
  best_r2 = subset(r2, m_score < m_score_cutoff)

  combined = merge(best_r1, best_r2, by=unique_merge_id)
  print(paste("After merging, there are", dim(combined)[1], "points left"))
  # remove lower peak groups and NAs
  df = subset(combined, peak_group_rank.x == 1 & peak_group_rank.y == 1)
  df = subset(df, ! is.na(RT.x) & ! is.na(RT.y) )
  df = subset(df, decoy.x == FALSE & decoy.y == FALSE)
  print(paste("After filtering, there are", dim(df)[1], "points left"))

  # smoothing spline model
  # -- the first value is the (expected) input and the second value is the
  #    (expected output). Since we want to predict how to convert from slave to
  #    master, slave is first and master is second
  sm = smooth.spline(df$RT.y,df$RT.x,cv=T)

  # print graph
  svg(filename)
  x_t = df$RT.y
  y_t = df$RT.x
  plot(x_t,y_t, cex=0.1, main="Plot identified peak groups", xlab="Run 1", ylab="Run 2")
  lines(sm, col="red", lwd=0.8)
  dev.off()
  return(sm)
}

# trailingOnly=TRUE means that only arguments after --args are returned
args <- commandArgs(trailingOnly = TRUE)

#
# Step 1: parse arguments, read input
# 
in_master <- args[1]
in_second <- args[2]
fdr_cutoff <- as.double(args[3])
rm(args)

if (in_master == in_second) 
{
  print("The two files are identical, will exit...")
  quit()
}
print("Will read input files")
r1 = read.csv(in_master,header=TRUE, sep="\t")
r2 = read.csv(in_second,header=TRUE, sep="\t")

#
# Step 2: calculate the spline
# 
print("Read input files...will start constructing the spline")
graph_filename = paste(strsplit(in_second, ".xls"), ".svg")
smooth.model = get_smoothing_spline(r1, r2, fdr_cutoff, graph_filename)

r1$aligned_rt = r1$RT
r2$aligned_rt = predict(smooth.model,r2$RT)$y 
r1$aligned_rt_diff = r1$RT - r1$aligned_rt
r2$aligned_rt_diff = r2$RT - r2$aligned_rt

#
# Step 3: Write output
# 
print("Constructed the spline, will write output files.")
write.table(r1, paste(strsplit(in_master, ".xls" ), "_aligned.xls", sep=""),sep="\t", row.names=FALSE)
write.table(r2, paste(strsplit(in_second, ".xls" ), "_aligned.xls", sep=""),sep="\t", row.names=FALSE)
print("Done, program finished successfully")

