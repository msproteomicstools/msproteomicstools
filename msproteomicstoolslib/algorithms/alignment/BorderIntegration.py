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
import msproteomicstoolslib.algorithms.graphs.graphs as graphs

def integrationBorderShortestPath(selected_pg, target_run, transformation_collection_, tree):
    """Determine the optimal integration border by using the shortest path in the MST

    Args:
        selected_pg(list(GeneralPeakGroup)): list of selected peakgroups (e.g. those passing the quality threshold)
        target_run(String): run id of the target run (where value is missing)
        transformation_collection_(:class:`.LightTransformationData`): structure to hold binary transformations between two different retention time spaces
        tree(list(tuple)): a minimum spanning tree (MST) represented as list of edges (for example [('0', '1'), ('1', '2')] ). Node names need to correspond to run ids.

    Returns:
        A tuple of (left_integration_border, right_integration_border)
    """

    # Get the shortest path through the MST
    available_runs = [pg.peptide.run.get_id() for pg in selected_pg if pg.get_fdr_score() < 1.0]
    final_path = graphs.findShortestMSTPath(tree, target_run, available_runs)

    # Extract the starting point and the starting borders (left/right)
    source_run = final_path[-1]
    best_pg = [pg for pg in selected_pg if pg.peptide.run.get_id() == source_run][0]
    rwidth = float(best_pg.get_value("rightWidth"))
    lwidth = float(best_pg.get_value("leftWidth"))
    # if verb: print "Compare from", source_run, "directly to", target_run,\
    #     transformation_collection_.getTrafo(source_run, target_run).predict([lwidth])[0], \
    #     transformation_collection_.getTrafo(source_run, target_run).predict([rwidth])[0]

    # Traverse path in reverse order
    for target_run_ in reversed(final_path[:-1]):

        # Transform, go to next step
        lwidth = transformation_collection_.getTrafo(source_run, target_run_).predict([lwidth])[0]
        rwidth = transformation_collection_.getTrafo(source_run, target_run_).predict([rwidth])[0]
        source_run = target_run_

    return lwidth, rwidth

def integrationBorderShortestDistance(selected_pg, target_run, transformation_collection_, mat, rmap):
    """Determine the optimal integration border by using the shortest distance (direct transformation)

    Args:
        selected_pg(list(GeneralPeakGroup)): list of selected peakgroups (e.g. those passing the quality threshold)
        target_run(String): run id of the target run (where value is missing)
        transformation_collection_(:class:`.LightTransformationData`): structure to hold binary transformations between two different retention time spaces
        mat(numpy matrix): distance matrix for all runs (as returned by algorithms.alignment.AlignmentMST.getDistanceMatrix)
        rmap(dict): mapping run ids to matrix columns (e.g. {"run_0" : 0, "run_1" : 1})

    Returns:
        A tuple of (left_integration_border, right_integration_border)
    """

    # Available runs which have a reliable RT value (no noise integration values ... )
    available_runs = [pg.peptide.run.get_id() for pg in selected_pg if pg.get_fdr_score() < 1.0]

    # Select curaent row from the matrix, reduce to columns of runs for which
    # we actually have a value and select closest among these
    current_matrix_row = mat[rmap[target_run],]
    current_matrix_row = [ (current_matrix_row[ rmap[curr] ], curr) for curr in available_runs]
    source_run = min(current_matrix_row)[1]

    # Transform
    best_pg = [pg for pg in selected_pg if pg.peptide.run.get_id() == source_run][0]
    rwidth = float(best_pg.get_value("rightWidth"))
    lwidth = float(best_pg.get_value("leftWidth"))
    leftW = transformation_collection_.getTrafo(source_run, target_run).predict([lwidth])[0]
    rightW = transformation_collection_.getTrafo(source_run, target_run).predict([rwidth])[0]
    return leftW, rightW

def integrationBorderReference(new_exp, selected_pg, rid, transformation_collection_, border_option):
    """Determine the optimal integration border by taking the mean of all other peakgroup boundaries using a reference run.

    Args:
        new_exp(AlignmentExperiment): experiment containing the aligned peakgroups
        selected_pg(list(GeneralPeakGroup)): list of selected peakgroups (e.g. those passing the quality threshold)
        rid(String): current run id
        transformation_collection_(:class:`.TransformationCollection`): specifying how to transform between retention times of different runs
        border_option(String): one of the following options ("mean", "median" "max_width"), determining how to aggregate multiple peak boundary information

    Returns:
        A tuple of (left_integration_border, right_integration_border) in the retention time space of the _reference_ run
    """
    current_run = [r for r in new_exp.runs if r.get_id() == rid][0]
    ref_id = transformation_collection_.getReferenceRunID()

    def convert_to_this(orig_runid, target_runid, ref_id, rt, transformation_collection_):
        """ Convert a retention time into one of the target RT space.
        
        Using the transformation collection
        """
        try:
            normalized_space_rt = transformation_collection_.getTransformation(orig_runid, ref_id).predict( [rt] )[0]
            return transformation_collection_.getTransformation(ref_id, target_runid).predict( [normalized_space_rt] )[0]
        except AttributeError as e:
            print("Could not convert from run %s to run %s (through reference run %s) -\
                    are you sure you gave the corresponding trafo file with \
                    the --in parameter?" % (orig_runid, target_runid, ref_id))
            print (e)
            raise e


    pg_lefts = []
    pg_rights = []
    for pg in selected_pg:

        rwidth = float(pg.get_value("rightWidth"))
        lwidth = float(pg.get_value("leftWidth"))
        this_run_rwidth = convert_to_this(pg.peptide.run.get_id(),
            current_run.get_id(), ref_id, rwidth, transformation_collection_)
        this_run_lwidth = convert_to_this(pg.peptide.run.get_id(),
            current_run.get_id(), ref_id, lwidth, transformation_collection_)

        pg_lefts.append(this_run_lwidth)
        pg_rights.append(this_run_rwidth)

    if border_option == "mean":
        integration_left = numpy.mean(pg_lefts)
        integration_right = numpy.mean(pg_rights)
    elif border_option == "median":
        integration_left = numpy.median(pg_lefts)
        integration_right = numpy.median(pg_rights)
    elif border_option == "max_width":
        integration_left = numpy.min(pg_lefts)
        integration_right = numpy.max(pg_rights)
    else:
        raise Exception("Unknown border determination option %s" % border_option)

    std_warning_level = 20
    if numpy.std(pg_rights) > std_warning_level or numpy.std(pg_lefts) > std_warning_level: 
        pass
        # TODO : what if they are not consistent ?
        # print "Std is too large", pg_rights, pg_lefts
        # print "overall mean, std", numpy.mean(pg_lefts), numpy.std(pg_lefts), numpy.mean(pg_rights), numpy.std(pg_rights)

    return integration_left, integration_right

