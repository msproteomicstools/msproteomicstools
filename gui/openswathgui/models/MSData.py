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
    from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader, inferMapping
    from msproteomicstoolslib.algorithms.alignment.MRExperiment import MRExperiment as Experiment
except ImportError:
    print "Could not find msproteomicstoolslib, certain functions are not available."

from ChromatogramTransition import ChromatogramTransition
from SwathRun import SwathRun
from SwathRunCollection import SwathRunCollection

REALIGN_RUNS = True
FDR_CUTOFF = 0.01
ONLY_SHOW_QUANTIFIED = False
ONLY_SHOW_QUANTIFIED = True

class PrecursorModel():
    """ A simplistic precursor model

    It is used internal to parse individual precursors. It can be initialized
    with an ID and knows how to parse the sequence and the charge from that
    string.

    Attributes:
        chrom_id(str): The chromatogram native id used to infer charge and peptide sequence
    """

    def __init__(self, chrom_id):
        self.chrom_id = chrom_id

    def getCharge(self):
        try:
            return int(self.chrom_id.split("/")[1].split("_")[0])
        except Exception:
            return "NA"

    def getFullSequence(self):
        try:
            assert self.chrom_id.find("/") != -1
            seq = self.chrom_id.split("/")[0].split("_")[-1]
            return seq
        except Exception:
            return self.chrom_id

class DataModel(object):
    """The main data model, provides access to all raw data

    It stores the references to individual :class:`.SwathRun` objects and can
    be initialized from a list of files. Each "load" method is
    responsible for setting the self.runs parameter.

    Attributes:
        runs(list of :class:`.SwathRun`): The MS runs which are handled by this class
        runs(bool): Whether to draw individual transitions
    """

    def __init__(self, fdr_cutoff = FDR_CUTOFF, only_quantified = ONLY_SHOW_QUANTIFIED):
        self.runs = []
        self.fdr_cutoff = fdr_cutoff
        self.only_show_quantified = only_quantified
        self.draw_transitions_ = False

    #
    ## Getters
    #
    def getStatus(self):
        """
        Returns
        -----------
        str:
            Returns its own status (number of transitions etc.) for the status bar.
        """
        if len(self.get_runs()) == 0:
            return "Ready"

        tr_cnt = 0
        for r in self.get_runs():
            tr_cnt += r.getTransitionCount()

        return '%s Transitions, %s Peptides (total %s Transitions)' % ( 
            self.get_runs()[0].getTransitionCount(), len(self.get_runs()[0]._sequences_mapping), tr_cnt )

    def get_precursor_tree(self):
        """
        Returns the data models precursor tree structure

        Returns a list of :class:`.ChromatogramTransition` root elements (rows)
        to display in the left side tree view. Each element may contain nested
        :class:`.ChromatogramTransition` elements (tree elements).

        Returns
        -----------
        list of :class:`.ChromatogramTransition`:
            Root element(s) for the peptide tree
        """
        return self._build_tree()

    def get_runs(self):
        """
        Returns the list of :class:`.SwathRun` objects of this current data model

        Returns
        -----------
        list of :class:`.SwathRun`
            The main content of the class is returned, its list of :class:`.SwathRun`
        """
        return self.runs

    def getDrawTransitions(self):
        """
        Returns
        -----------
        bool
            Whether to draw transitions or not
        """
        return self.draw_transitions_

    def setDrawTransitions(self, draw_transitions):
        """
        Whether to draw individual transitions or not
        """
        self.draw_transitions_ = draw_transitions

    #
    ## Data loading
    #
    def loadSqMassFiles(self, filenames):

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_sql(filenames)
        self.runs = [run for run in swathfiles.getSwathFiles()]

    def loadFiles(self, filenames):
        """
        Load a set of chromatogram files (no peakgroup information).

        Args:
            filenames(list of str): List of filepaths containing the chromatograms
        """

        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_files(filenames)
        self.runs = [run for run in swathfiles.getSwathFiles()]

    def loadMixedFiles(self, rawdata_files, aligned_pg_files, fileType):
        """ Load files that contain raw data files and aligned peakgroup files.

        Since no mapping is present here, we need to infer it from the data.
        Basically, we try to map the column align_runid to the filenames of the
        input .chrom.mzML hoping that the user did not change the filenames.

        Parameters
        ----------
        rawdata_files : list of str
            List of paths to chrom.mzML files
        aligned_pg_files : list of str
            List of paths to output files of the FeatureAligner
        fileType : str
            Description of the type of file the metadata file (valid: simple, traml, openswath)
        """

        print "Input contained no mapping of run_id to the chromatograms."
        print "Try to infer mapping for filetype %s - if this fails, please provide a yaml input." % fileType

        precursors_mapping = {}
        sequences_mapping = {}
        protein_mapping = {}
        mapping = {}
        inferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping, sequences_mapping, protein_mapping, fileType=fileType)
        print "Found the following mapping: mapping", mapping

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        if fileType == "sqmass":
            swathfiles.initialize_from_sql_map(mapping, rawdata_files, precursors_mapping, sequences_mapping, protein_mapping)
        elif self.only_show_quantified:
            swathfiles.initialize_from_chromatograms(mapping, precursors_mapping, sequences_mapping, protein_mapping)
        else:
            swathfiles.initialize_from_chromatograms(mapping)
        self.runs = [run for run in swathfiles.getSwathFiles()]

        if not fileType in ["simple", "traml"]:
            self._read_peakgroup_files(aligned_pg_files, swathfiles)

        print "Find in total a collection of %s runs." % len(swathfiles.getRunIds() )
                    
    def load_from_yaml(self, yamlfile):
        """
        Load a yaml file containing a mapping of chromatogram files and aligned peakgroup files.

        Parameters
        ----------
        yamlfile : str
            Filepath to the yaml file for loading
        """

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
            - assay_rt
            - Intensity
            - align_runid
            - transition_group_id
        """

        # Read in the peakgroup files, parse them and map across runs
        reader = SWATHScoringReader.newReader(aligned_pg_files, "openswath", readmethod="gui", errorHandling="loose")
        new_exp = Experiment()
        new_exp.runs = reader.parse_files(REALIGN_RUNS)
        multipeptides = new_exp.get_all_multipeptides(self.fdr_cutoff, verbose=False)

        # Build map of the PeptideName/Charge to the individual multipeptide
        peakgroup_map = {}
        for m in multipeptides:
            pg = m.find_best_peptide_pg()
            identifier = pg.get_value("FullPeptideName") + "/" + pg.get_value("Charge")
            # identifier for precursor, see msproteomicstoolslib/format/SWATHScoringMapper.py
            peakgroup_map[ identifier ] = m
            peakgroup_map[ identifier + "_pr" ] = m

        for swathrun in swathfiles.getSwathFiles():
            if self.only_show_quantified:
                intersection = set(swathrun.get_all_precursor_ids()).intersection( peakgroup_map.keys() )
                todelete = set(swathrun.get_all_precursor_ids()).difference(intersection)
                if len(intersection) == 0:
                    print "Could not find any intersection between identifiers in your transition file and the provided chromatograms"
                    print len(intersection)
                swathrun.remove_precursors(todelete)

            # for each precursor in this run, identify the best peakgroup and store the value
            for precursor_id in swathrun.get_all_precursor_ids():
                if not peakgroup_map.has_key(precursor_id): 
                    continue

                m = peakgroup_map[ precursor_id ]
                if m.hasPrecursorGroup(swathrun.runid):
                    for pg in m.getPrecursorGroup(swathrun.runid).getAllPeakgroups():
                        l,r       = [ float(pg.get_value("leftWidth")), float(pg.get_value("rightWidth")) ]
                        fdrscore  = float(pg.get_value("m_score"))
                        assay_rt  = float(pg.get_value("assay_rt"))
                        intensity = float(pg.get_value("Intensity"))
                        swathrun.add_peakgroup_data(precursor_id, l, r, fdrscore, intensity, assay_rt)

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
        """
        Build tree of :class:`.ChromatogramTransition` objects for display (see get_precursor_tree)
        """

        peptide_sequences = set([])
        for r in self.get_runs():
            peptide_sequences.update( r.get_all_peptide_sequences() )

        proteins = set([])
        for r in self.get_runs():
            proteins.update( r.get_all_proteins() )

        if len(proteins) == 0:
            elements = self._build_tree_sequences(peptide_sequences)
        else:

            elements = []
            for protein in proteins:

                # get all sequences from all runs
                sequences = set([])
                for r in self.get_runs():
                    sequences.update( r.get_sequence_for_protein(protein) )

                pelements = self._build_tree_sequences(sequences)

                elements.append(ChromatogramTransition(protein,
                                                       "NA",
                                                       pelements, 
                                                       datatype="Protein") )

        return elements

    def _build_tree_sequences(self, peptide_sequences):
        """
        Build tree of :class:`.ChromatogramTransition` objects for display (see get_precursor_tree)
        """
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

                    # Only add transition data if individual transitions should be drawn
                    if self.draw_transitions_:
                        tr_elements.append(
                            ChromatogramTransition(tr,
                                                   -1,
                                                   [],
                                                   fullName=tr,
                                                   peptideSequence = pm.getFullSequence(),
                                                   datatype="Transition") )

                pelements.append(ChromatogramTransition(p,
                                                        pm.getCharge(),
                                                        tr_elements, 
                                                        peptideSequence = pm.getFullSequence(),
                                                        datatype="Precursor") )

            elements.append(ChromatogramTransition(seq,
                                                   "NA",
                                                   pelements, 
                                                   datatype="Peptide",
                                                   peptideSequence=pm.getFullSequence()) )
        return elements

