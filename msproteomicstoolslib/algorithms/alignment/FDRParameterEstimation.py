#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013, ETH Zurich
For a full list of authors, refer to the file AUTHORS.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""

from __future__ import print_function
import numpy

class ParamEst(object):
    """
    Parameter estimation object

    In a first step the percentage of decoys of all peakgroups at the target
    fdr is computed (which is then taken as the "aim"). For this "aim" of decoy
    percentage, the class will try to estimate an fdr_cutoff such that the
    percentage of decoy precursors in the final reported result will correspond
    to the "aim". 

    If the parameter min_runs (at initialization) is higher than 1, only
    precursors that are identified in min_runs above the fdr_cutoff will be
    reported.

    >>> p = ParamEst()
    >>> decoy_frac = p.compute_decoy_frac(multipeptides, target_fdr)
    >>> fdr_cutoff_calculated = p.find_iterate_fdr(multipeptides, decoy_frac)
    """

    def __init__(self, min_runs=1, max_recursion=10, verbose=False):
        self.verbose = verbose
        self.min_runs = min_runs
        self.max_recursion = max_recursion

    def find_iterate_fdr(self, multipeptides, decoy_frac, recursion=0):
        """Iteratively find an q-value cutoff to reach the specified decoy fraction

        This function will step through multiple q-value thresholds and
        evaluate how many peptides have at least one peakgroup whose q-value
        (using get_fdr_score) is below that threshold. The q-value is then
        adapted until the fraction of decoys in the result is equal to the
        specified fraction given as input.

        Args:
            multipeptides(list(Multipeptide)): a list of
                multipeptides on which the alignment should be performed. After
                alignment, each peakgroup that should be quantified can be
                retrieved by calling get_selected_peakgroups() on the multipeptide.
            decoy_frac(float): fraction of decoys that should be present in the result (e.g. 0.01)
            recursion(int): recursion level (should be set to zero)

        Returns:
            None
        """


        # Starting guess
        start = 0.05 / (10**recursion)
        end = 1.0 / (10**recursion)
        stepsize = start

        if self.verbose: 
            print("Recurse", recursion)

        if recursion > self.max_recursion:
            raise Exception("Recursed too much in FDR iteration.")

        decoy_pcnt = decoy_frac*100
        val_005 = self._calc_precursor_fr(multipeptides, (start+stepsize)/100.0 )*100
        val_1 = self._calc_precursor_fr(multipeptides, end/100.0 )*100

        if self.verbose: print("Decoy pcnt aim:", decoy_pcnt)
        if self.verbose: print("Aim, high_value, low_value", decoy_pcnt, val_1, val_005)

        # Check if 
        # i)   we already found the correct value (no need for iteration)
        # ii)  our computed value lies between 0.05% and 1% FDR cutoff
        # iii) is higher than 1% 
        if recursion == 0 and abs(decoy_pcnt - val_1) < 1e-6:
            return decoy_frac

        elif decoy_pcnt < val_005:
            return self.find_iterate_fdr(multipeptides, decoy_frac, recursion=recursion+1)

        elif decoy_pcnt > val_1:
            if self.verbose: 
                print("choose larger step from 0.5 on")

            start = 0.5
            end = 100.0
            stepsize = 0.5
            if recursion > 1: # pragma: no cover
                raise Exception ("Decreased start / end but the values was too large? Should never happen.")
        else:
            # All is fine, we are within the limits
            pass

        fdrrange = numpy.arange(start, end + 2*stepsize, stepsize) # add 2 extra steps for edge cases
        return self._find_iterate_fdr(multipeptides, decoy_frac, fdrrange)

    def _find_iterate_fdr(self, multipeptides, decoy_frac, fdrrange):
        """
        Iterate through the given range of q-value cutoffs.
        
        Stop if the calculated decoy fraction is larger than the decoy_fraction
        parameter.
        """
        decoy_pcnt = decoy_frac*100

        if self.verbose: 
            print("mScore_cutoff", "Calc-precursor-FDR")

        for fdr in fdrrange:
            calc_fdr = self._calc_precursor_fr(multipeptides, fdr/100.0 )*100
            if self.verbose: 
                print(fdr, calc_fdr)

            # Break if the calculated FDR is higher than our aim -> the true
            # value lies between this and the previous value.
            if calc_fdr > decoy_pcnt:
                break

            # Break if the calculated FDR is (almost) exactly our aim -> the
            # user wants to have exactly this FDR
            if abs(calc_fdr - decoy_pcnt) < 1e-6:
                break

            prev_fdr = fdr
            prev_calc_fdr = calc_fdr

        # The last value is extremely close to the true one
        if abs(calc_fdr - decoy_pcnt) < 1e-6:
            return fdr/100.0

        # We have run through without stopping
        if prev_fdr == fdr:
            raise Exception("Parameter estimation did not reach a high enough value")

        # Linear interpolation
        res = prev_fdr + (fdr-prev_fdr) * (decoy_pcnt-prev_calc_fdr)/(calc_fdr-prev_calc_fdr)
        return res/100.0

    def _calc_precursor_fr(self, multipeptides, target_fdr):
        """ Calculate how many of the *precursors* are decoy for a given cutoff.
        """
        min_runs = self.min_runs
        allpg_cnt = 0
        alldecoypg_cnt = 0
        for mpep in multipeptides:
            count = 0
            decoy = False
            for prgr in mpep.getPrecursorGroups():
                for pep in prgr:
                  if pep.get_best_peakgroup().get_fdr_score() < target_fdr:
                      count += 1
                  if pep.get_decoy():
                      decoy = True
            if count >= min_runs:
                allpg_cnt += 1
            if decoy and count >= min_runs:
                alldecoypg_cnt += 1

        # Check for degenerate case
        if allpg_cnt == 0:
            return 0.0

        return alldecoypg_cnt *1.0 / allpg_cnt

    def compute_decoy_frac(self, multipeptides, target_fdr):
        """ Calculate how many of the *peakgroups* are decoy for a given cutoff.

        This provides an estimator of the ration of false positives to decoy
        peak groups which the mProphet / pyProphet algorithm has previously
        estimated using the Storey-Tibshirani method.
        """
        allpg_cnt = 0
        alldecoypg_cnt = 0
        for mpep in multipeptides:
            for prgr in mpep.getPrecursorGroups():
                for pep in prgr:
                    if pep.get_best_peakgroup().get_fdr_score() < target_fdr:
                        allpg_cnt += 1
                        if pep.get_decoy():
                            alldecoypg_cnt += 1

        decoy_frac = alldecoypg_cnt *1.0 / allpg_cnt
        return decoy_frac

