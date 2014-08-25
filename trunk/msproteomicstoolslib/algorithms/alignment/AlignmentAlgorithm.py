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

from sys import stdout
from msproteomicstoolslib.algorithms.alignment.Multipeptide import Multipeptide
import numpy

class Cluster:
    """ A helper class representation of a cluster (used in AlignmentAlgorithm)
    """

    def __init__(self, peakgroups):
        """
        Initialization 

        Args:
            peakgroups(list(SWATHScoringReader.MinimalPeakGroup))
        """
        self.peakgroups = peakgroups

    def select_one_per_run(self, verbose=False):
      """Ensure that only one peakgroup is selected per run.

      If there are multiple peakgroups selected, only the best one is retained.
      """
      run_ids = {}
      if verbose: print "len pg ", len(self.peakgroups)
      for pg in self.peakgroups:
          rid = pg.peptide.get_run_id()
          if rid in run_ids:
              if verbose: print "have run id", rid, "multiple times", pg.get_fdr_score(),\
                 "/", pg.get_normalized_retentiontime(), " vs ", \
                 run_ids[rid].get_fdr_score(), "/", run_ids[rid].get_normalized_retentiontime()
              if run_ids[rid].get_fdr_score() > pg.get_fdr_score():
                  run_ids[rid] = pg
          else: run_ids[rid] = pg
      self.peakgroups = run_ids.values()

    def get_total_score(self):
      """
      Calculate the total score of a cluster (multiplication of probabilities)
      """
      mult = 1
      for pg in self.peakgroups:
        mult = mult * pg.get_fdr_score()
      return mult

    def get_median_rt(self):
      """
      Calculate the total score of a cluster (multiplication of probabilities)
      """

      return numpy.median( [pg.get_normalized_retentiontime() for pg in self.peakgroups] )

class AlignmentAlgorithm():
    """ A class of alignment algorithms
    """
    def __init__(self): 
        self.verbose = False

    def align_features(self, multipeptides, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, method="best_overall"):
        """ Perform the alignment on a set of multipeptides

        Args:
            multipeptides(list(Multipeptide)): a list of
                multipeptides on which the alignment should be performed. After
                alignment, each peakgroup that should be quantified can be
                retrieved by calling get_selected_peakgroups() on the multipeptide.
            rt_diff_cutoff(float): maximal allowed RT difference used in the clustering or the 
            fdr_cutoff(float): FDR cutoff for "seeds" (e.g. the "known good" cutoff, for example 0.01)
            aligned_fdr_cutoff(float): maximal FDR cutoff to still be considered for the clustering
            method(String): either best_overall or best_cluster_score (or global_best_cluster_score, global_best_overall)

        Returns:
            Alignment object with alignment statistics
        """

        class AlignmentStatistics():
            def __init__(self): 
                self.nr_aligned = 0
                self.nr_changed = 0
                self.nr_quantified = 0
                self.could_not_align = 0
                self.fdr_good_in_all = 0
        astas = AlignmentStatistics()

        for mpep in multipeptides:

            if mpep.all_above_cutoff(fdr_cutoff) and method in ["best_cluster_score", "best_overall"]:
                # In non-global algorithms, do not re-align if the fdr is above the threshold in all runs
                for p in mpep.get_peptides():
                    p.get_best_peakgroup().select_this_peakgroup()
                    astas.nr_quantified += 1
                continue

            if method == "global_best_cluster_score" or method == "best_cluster_score":
                self.align_features_cluster_(mpep, astas, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, method)
            elif method == "global_best_overall" or method == "best_overall":
                self.align_features_best_(mpep, astas, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, method)
            else:
                raise Exception("Method '%s' unknown" % method)
        return astas

    def align_features_cluster_(self, m, a, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, method):
        verb=self.verbose
        if verb: print "00000000000000000000000000000000000 new peptide ", m.get_peptides()[0].sequence

        # i) get all RTs above the cutoff
        for p in m.get_peptides(): # loop over runs
            pg = p.get_best_peakgroup()
            if verb: print "best rt", pg.get_normalized_retentiontime(), pg.peptide.run.get_id(), pg.get_fdr_score()
        
        groups = [ pg 
            for p in m.get_peptides() # loop over all peptides
                for pg in p.get_all_peakgroups() # loop over all peakgroups
                    if pg.get_fdr_score() < aligned_fdr_cutoff
        ]

        # Check for empty groups
        if len(groups) == 0:
            return

        # do the clustering
        from cluster import HierarchicalClustering
        cl = HierarchicalClustering(groups, lambda x,y: abs(x.get_normalized_retentiontime()-y.get_normalized_retentiontime()))
        clusters_rt = cl.getlevel(rt_diff_cutoff) # for large clusters, this is the the bottleneck! 
        clusters_rt_obj = [Cluster(c) for c in clusters_rt]
        # if there was only one group, we need to prepare a special object of size one
        if len(groups) == 1: clusters_rt_obj = [Cluster( groups )]

        if verb: print "==== Clusters "
        # make sure only one is selected from each run...
        for c in clusters_rt_obj: 
            c.select_one_per_run(self.verbose)
            if verb:
                print " - Cluster with score", c.get_total_score(), "at", c.get_median_rt()
                for pg in c.peakgroups: print pg.get_normalized_retentiontime(), pg.peptide.run.get_id()
          
        if len(clusters_rt_obj) == 1 :
            # great, only one cluster
            bestcluster = clusters_rt_obj[0]
        else:
            # select the best cluster
            bestcluster = clusters_rt_obj[0] # self.getBestCluster(clusters_rt, rt_maximal_distance, method = clusterMode)
            for cluster in clusters_rt_obj:
                if cluster.get_total_score() < bestcluster.get_total_score(): bestcluster = cluster

        for pg in bestcluster.peakgroups:
          a.nr_quantified += 1
          pg.select_this_peakgroup()
          if pg.get_fdr_score() > fdr_cutoff:
              a.nr_aligned += 1
              if pg.get_normalized_retentiontime() != pg.peptide.get_best_peakgroup().get_normalized_retentiontime():
                  a.nr_changed += 1
                  if verb: print "FDR new align", pg.peptide.get_best_peakgroup().print_out(), "\tnew ====> ", pg.print_out()
              else:
                  if verb: print "FDR boost", pg.peptide.get_best_peakgroup().print_out(), " old ====> ", pg.print_out()
          else:
            if verb: print "no need to align", pg.print_out()

    def align_features_best_(self, m, a, rt_diff_cutoff, fdr_cutoff, aligned_fdr_cutoff, method):
        verb=self.verbose
        if verb: print "00000000000000000000000000000000000 new peptide (best overall)", m.get_peptides()[0].sequence

        # If we just choose the cluster with the "best" peptide, we find find the best peptide over all runs
        best = m.find_best_peptide_pg()
        best_rt_diff = best.get_normalized_retentiontime()
        if verb: print "=====\nFDR best", best.print_out()

        for p in m.get_peptides(): # loop over runs
            current_best_pg = p.get_best_peakgroup()

            if current_best_pg.get_fdr_score() < fdr_cutoff and method == "best_overall":
                # The pg is below the fdr cutoff, we just take it and go with it in "best overall"
                current_best_pg.select_this_peakgroup()
                a.nr_quantified += 1
                if verb: 
                    print " == Selected peakgroup ", current_best_pg.print_out()

                continue

            # #################################################################
            # In this run, the peptide is above the FDR cutoff. We will now:
            #   - find all peakgroups that are within the retention time cutoff (continue if none are found)
            #   - of those select the peakgroup with the best score
            #   - if the best-scoring peakgroup is acceptable (<aligned_fdr_cutoff), mark it as selected

            matching_peakgroups = [pg_ for pg_ in p.get_all_peakgroups() 
                if (abs(float(pg_.get_normalized_retentiontime()) - float(best_rt_diff)) < rt_diff_cutoff)]

            if len(matching_peakgroups) == 0: 
                continue

            bestScoringPG = min(matching_peakgroups, key=lambda x: float(x.get_fdr_score()))

            if bestScoringPG.get_fdr_score() < aligned_fdr_cutoff: 
                bestScoringPG.select_this_peakgroup()

                a.nr_aligned += 1
                a.nr_quantified += 1
                if current_best_pg.get_normalized_retentiontime() != bestScoringPG.get_normalized_retentiontime():
                  a.nr_changed += 1
                  if verb: print "FDR new align", current_best_pg.print_out(), "\tnew ====> ", bestScoringPG.print_out()
                else:
                  if verb: print "FDR boost", current_best_pg.print_out(), " old ====> ", bestScoringPG.print_out()
            else:
                if verb: print "could not align", current_best_pg.peptide.run.get_id(), current_best_pg.peptide.run.orig_filename, "best rt_diff was ", \
                      abs(float(bestScoringPG.get_normalized_retentiontime()) - float(best_rt_diff)), "best score", \
                      bestScoringPG.get_fdr_score() 
                a.could_not_align += 1

