#!/usr/bin/python
# -*- coding: utf-8 -*-
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

import os

try:
    from msproteomicstoolslib.format.TransformationCollection import TransformationCollection
    from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader
    from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment as Experiment
except ImportError:
    print "Could not find msproteomicstoolslib, certain functions are not available."

from ChromatogramTransition import ChromatogramTransition
from SwathRun import SwathRun
from SwathRunCollection import SwathRunCollection

DRAW_TRANSITIONS = False
REALIGN_RUNS = True
FDR_CUTOFF = 0.01
ONLY_SHOW_QUANTIFIED = False
ONLY_SHOW_QUANTIFIED = True

class PrecursorModel():
    """ A simplistic precursor model

    It is initialized with an ID and knows how to parse the sequence and the
    charge from that string.
    """

    def __init__(self, chrom_id):
        self.chrom_id = chrom_id

    def getCharge(self):
        try:
            return self.chrom_id.split("/")[1].split("_")[0]
        except Exception:
            return "NA"

    def getFullSequence(self):
        try:
            return self.chrom_id.split("/")[0].split("_")[-1]
        except Exception:
            return "NA"

class DataModel(object):
    """The main data model

    It stores the SwathRun objects and can be initialized from a list of files. 

    Attributes:
        runs:     A list of SwathRun objects

    Methods:
        get_precursors_tree:  Returns a list of ChromatogramTransition root elements (rows) to display in the left side tree
                              view. Each element may contain nested ChromatogramTransition elements (tree elements).

        get_runs:             Returns the list of SwathRun objects of this
                              current data model

        loadFiles:            Load a set of chromatogram files (no peakgroup information).

        loadMixedFiles:       Load a set of chromatogram files (flat file list) and aligned peakgroup files. 

        load_from_yaml:       Load a yaml file containing a mapping of chromatogram files and aligned peakgroup files.

    """

    def __init__(self):
        self.runs = []

    #
    ## Getters
    #
    def getStatus(self):
        if len(self.get_runs()) == 0:
            return "Ready"

        tr_cnt = 0
        for r in self.get_runs():
            tr_cnt += r.getTransitionCount()

        return '%s Transitions, %s Peptides (total %s Transitions)' % ( 
            self.get_runs()[0].getTransitionCount(), len(self.get_runs()[0]._sequences_mapping), tr_cnt )

    def get_precursor_tree(self):
        return self._build_tree()

    def get_runs(self):
        return self.runs

    #
    ## Data loading
    #
    def loadFiles(self, filenames):

        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_files(filenames)
        self.runs = [run for run in swathfiles.getSwathFiles()]

    def loadMixedFiles(self, rawdata_files, aligned_pg_files):
        """ Load files that contain raw data files and aligned peakgroup files.

        Since no mapping is present here, we need to infer it from the data.
        Basically, we try to map the column align_runid to the filenames of the
        input .chrom.mzML hoping that the user did not change the filenames.
        """

        print "Input contained no mapping of run_id to the chromatograms."
        print "Try to infer mapping - if this fails, please provide a yaml input."
        mapping = {}
        import csv, os
        for file_nr, f in enumerate(aligned_pg_files):
            header_dict = {}
            if f.endswith('.gz'):
                import gzip 
                filehandler = gzip.open(f,'rb')
            else:
                filehandler = open(f)
            reader = csv.reader(filehandler, delimiter="\t")
            header = reader.next()
            for i,n in enumerate(header):
                header_dict[n] = i
            if not header_dict.has_key("align_origfilename") or not header_dict.has_key("align_runid"):
                print header_dict
                raise Exception("need column header align_origfilename and align_runid")

            for this_row in reader:

                # 1. Get the original filename (find a non-NA entry) and the corresponding run id
                if len(this_row) == 0: continue
                if this_row[ header_dict["align_origfilename"] ] == "NA": continue
                aligned_id = os.path.basename(this_row[ header_dict["align_runid"] ])
                if aligned_id in mapping: continue 
                aligned_fname = os.path.basename(this_row[ header_dict["align_origfilename"] ])

                # 2. Go through all chromatogram input files and try to find
                # one that matches the one from align_origfilename
                for rfile in rawdata_files:
                    # 2.1 remove common file endings from the raw data
                    rfile_base = os.path.basename(rfile)
                    for ending in [".mzML", ".chrom"]:
                        rfile_base = rfile_base.split(ending)[0]
                    # 2.2 remove common file endings from the tsv data
                    for ending in [".tsv", ".csv", ".xls", "_with_dscore", "_all_peakgroups"]:
                        aligned_fname = aligned_fname.split(ending)[0]
                    # 2.3 Check if we have a match
                    if aligned_fname == rfile_base:
                        print "- Found match:", rfile_base, aligned_fname
                        mapping[aligned_id] = [rfile]

        print "Found the following mapping: mapping", mapping

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_chromatograms(mapping)
        self.runs = [run for run in swathfiles.getSwathFiles()]
        self._read_peakgroup_files(aligned_pg_files, swathfiles)
        print "Find in total a collection of %s runs." % len(swathfiles.getRunIds() )
                    
    def load_from_yaml(self, yamlfile):

        import yaml
        data = yaml.load(open(yamlfile) )["AlignedSwathRuns"]
        alignment_files = data["PeakGroupData"]
        trafo_fnames = [d["trafo_file"] for d in data["RawData"]]
        self._loadFiles_with_peakgroups(data["RawData"], alignment_files)

    #
    ## Private functions
    #
    def _read_peakgroup_files(self, aligned_pg_files, swathfiles):
        """
        The peakgroup files have to have the following columns:
            - FullPeptideName
            - Charge
            - leftWidth
            - rightWidth
            - m_score
            - Intensity
            - align_runid
            - transition_group_id
        """

        # Read in the peakgroup files, parse them and map across runs
        reader = SWATHScoringReader.newReader(aligned_pg_files, "openswath", readmethod="gui", errorHandling="loose")
        new_exp = Experiment()
        new_exp.runs = reader.parse_files(REALIGN_RUNS)
        multipeptides = new_exp.get_all_multipeptides(FDR_CUTOFF, verbose=False)

        # Build map of the PeptideName/Charge to the individual multipeptide
        peakgroup_map = {}
        for m in multipeptides:
            pg = m.find_best_peptide_pg()
            identifier = pg.get_value("FullPeptideName") + "/" + pg.get_value("Charge")
            peakgroup_map[ identifier ] = m

        for swathrun in swathfiles.getSwathFiles():
            if ONLY_SHOW_QUANTIFIED:
                intersection = set(swathrun.get_all_precursor_ids()).intersection( peakgroup_map.keys() )
                todelete = set(swathrun.get_all_precursor_ids()).difference(intersection)
                swathrun.remove_precursors(todelete)

            # for each precursor in this run, identify the best peakgroup and store the value
            for precursor_id in swathrun.get_all_precursor_ids():
                if not peakgroup_map.has_key(precursor_id): 
                    continue

                m = peakgroup_map[ precursor_id ]
                if m.has_peptide(swathrun.runid):
                    pg = m.get_peptide(swathrun.runid).get_best_peakgroup()
                    swathrun._range_mapping[precursor_id]       = [ float(pg.get_value("leftWidth")), float(pg.get_value("rightWidth")) ]
                    swathrun._score_mapping[precursor_id]       = float(pg.get_value("m_score"))
                    swathrun._intensity_mapping[precursor_id]   = float(pg.get_value("Intensity"))

    def _loadFiles_with_peakgroups(self, RawData, aligned_pg_files):

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        try:
            swathfiles.initialize_from_directories( dict( [ (d["id"], d["directory"]) for d in RawData] ) )
        except KeyError:
            swathfiles.initialize_from_chromatograms( dict( [ (d["id"], d["chromatograms"]) for d in RawData] ) )
        self.runs = [run for run in swathfiles.getSwathFiles()]
        print "Find in total a collection of %s runs." % len(swathfiles.getRunIds() )

        try:
            self._read_trafo(RawData)
        except IOError:
            self._read_peakgroup_files(aligned_pg_files, swathfiles)

    def _read_trafo(self, trafo_filenames):
        # Read the transformations
        transformation_collection_ = TransformationCollection()
        for filename in [d["trafo_file"] for d in trafo_filenames]:
          transformation_collection_.readTransformationData(filename)
        transformation_collection_.initialize_from_data(reverse=True)

    def _build_tree(self):

        peptide_sequences = set([])
        for r in self.get_runs():
            peptide_sequences.update( r.get_all_peptide_sequences() )

        ## The peptide sequences are our top-level items
        # print "pepseqs", peptide_sequences
        elements = []
        for seq in peptide_sequences:

            # get all precursors from all runs
            precursors = set([])
            for r in self.get_runs():
                precursors.update( r.get_precursors_for_sequence(seq) )

            # print "found precursros", precursors
            pelements = []
            for p in precursors:

                # get all transitions from all runs
                transitions = set([])
                for r in self.get_runs():
                    transitions.update( r.get_transitions_for_precursor(p) )
                tr_elements = []
                pm = PrecursorModel(p)
                for tr in transitions:
                    if DRAW_TRANSITIONS:
                        tr_elements.append(ChromatogramTransition(tr, -1, [], fullName=tr, peptideSequence = pm.getFullSequence(), datatype="Transition") )
                pelements.append(ChromatogramTransition(p, pm.getCharge(), tr_elements, 
                       peptideSequence = pm.getFullSequence(), datatype="Precursor") )
            elements.append(ChromatogramTransition(seq, "NA", pelements, datatype="Peptide", 
                       peptideSequence=pm.getFullSequence()) )
        return elements

