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

DRAW_TRANSITIONS = False
REALIGN_RUNS = True
FDR_CUTOFF = 0.01
ONLY_SHOW_QUANTIFIED = False
ONLY_SHOW_QUANTIFIED = True

class RunDataModel():
    """Data Model for a single file from one run
    """

    def __init__(self, run, filename, load_in_memory=False):
        import os
        self._run = run
        self._filename = filename
        self._basename = os.path.basename(filename)

        # Map which holds the relationship between a precursor and its
        # corresponding transitions.
        self._precursor_mapping = {}
        # Map which holds the relationship between a sequence and its precursors
        self._sequences_mapping = {}

        # may contain extra info for each precursor
        self._range_mapping = {}
        self._score_mapping = {}
        self._intensity_mapping = {}

        if load_in_memory:
            # Load 10s, 580 MB for 3 files. takes ca 0.15 seconds to display
            self._group_by_precursor_in_memory()
            self._in_memory = True
        else:
            # Load 0.7 s, 55.9 MB for 3 files. takes ca 0.15 seconds to display
            self._group_by_precursor()
            self._in_memory = False

        self._group_precursors_by_sequence()

    #
    ## Initialization
    #
    def _group_by_precursor(self):
        """
        Populate the mapping between precursors and the chromatogram ids.

        The precursor is of type 'PEPT[xx]IDE/3' with an optional (DECOY) tag
        in front. Different modifications will generate different precursors,
        but not different charge states.
        """
        
        openswath_format = self._has_openswath_format(self._run)
        if openswath_format:
            if len( self._run.info['offsets'] ) > 0:
                for key in self._run.info['offsets'].keys():
                    # specific to pymzl, we need to get rid of those two entries or
                    # when the key is zero
                    if key in ("indexList", "TIC"): continue
                    if len(key) == 0: continue
                    #
                    trgr_nr = self._compute_transitiongroup_from_key(key)
                    tmp = self._precursor_mapping.get(trgr_nr, [])
                    tmp.append(key)
                    self._precursor_mapping[trgr_nr] = tmp

        else:
            print "Fall-back: Could not group chromatograms by their id alone, try using precursor information."
            self._group_by_precursor_by_mass()

    def _group_by_precursor_by_mass(self):
        """
        Populate the mapping between precursors and the chromatogram ids using the chromatogram information.

        Try to use the mass or the precursor tag to infer which chromatograms
        belong to the same peptide. See _compute_transitiongroup_from_precursor
        """

        for key in self._run.info['offsets'].keys():
            # specific to pymzl, we need to get rid of those two entries or
            # when the key is zero
            if key in ("indexList", "TIC"): continue
            if len(key) == 0: continue
            #
            chromatogram = self._run[key]
            if len(chromatogram["precursors"]) == 1:
                precursor = chromatogram["precursors"][0]
                trgr_nr = self._compute_transitiongroup_from_precursor(precursor)
                tmp = self._precursor_mapping.get(trgr_nr, [])
                tmp.append(key)
                self._precursor_mapping[trgr_nr] = tmp

            else:
                print "Something is wrong"
                raise Exception("Could not parse chromatogram ids ... neither ids nor precursors lead to sensible grouping")

    def _group_by_precursor_in_memory(self):
        """
        Same as _group_by_precursor but loads all data in memory.
        """
        openswath_format = self._has_openswath_format(self._run)
        result = {}
        if openswath_format:
            for chromatogram in self._run:
                    key = chromatogram["id"]
                    trgr_nr = self._compute_transitiongroup_from_key(key)
                    tmp = self._precursor_mapping.get(trgr_nr, [])
                    tmp.append(key)
                    self._precursor_mapping[trgr_nr] = tmp

                    import copy
                    cc = copy.copy(chromatogram)
                    result[ key ] = cc
        else:
            raise Exception("Could not parse chromatogram ids ... ")
        
        # we work in memory
        self._run = result

    def _compute_transitiongroup_from_key(self, key):
        """ Transforms an input chromatogram id to a string of [DECOY]_xx/yy

        Possible Input Formats:

        i) [DECOY_]id_xx/yy_zz
        ii) [DECOY_]id_xx_yy

        Where the DECOY_ prefix is optional, xx is the sequence and yy is the
        charge. zz is an additional, optional annotation.

        We want to return only [DECOY]_xx/yy (ensuring that decoys get a
        different identifier than targets).
        """
        components = key.split("_")
        trgr_nr = str(components[1])
        if components[0].startswith("DECOY"):
            trgr_nr = "DECOY_" + str(components[2])

        # Format ii) (second component doesnt contain a slash)
        if trgr_nr.find("/") == -1:
            trgr_nr += "/" + str(components[-1])
        return trgr_nr

    def _compute_transitiongroup_from_precursor(self, precursor):
        """ Transforms an input precursor to a string of [DECOY]_xx/yy
        """
        charge = "0"
        if precursor.has_key("charge"):
            charge = str(precursor["charge"])

        # Check if there is a user-parameter describing the peptide sequence
        # which would allow us to asseble the SEQ/ch pair.
        if precursor.has_key("userParams"):
            if precursor["userParams"].has_key("peptide_sequence"):
                return precursor["userParams"]["peptide_sequence"] + "/" + charge
            if precursor["userParams"].has_key("PeptideSequence"):
                return precursor["userParams"]["PeptideSequence"] + "/" + charge

        if precursor.has_key("mz"):
            # TODO have some reproducible rounding here
            return str(precursor["mz"]) + "GENERIC/" + charge

    def _has_openswath_format(self, run):
        """Checks whether the chromatogram id follows a specific format which
        could allow to map chromatograms to precursors without reading the
        whole file.

        Namely, the format is expected to be [DECOY_]\d*_.*_.* from which one
        can infer that it is openswath format.
        """

        openswath_format = False
        if len( run.info['offsets'] ) > 0:
            keys = run.info['offsets'].keys()
            for key in run.info['offsets'].keys():
                if key in ("indexList", "TIC"): continue
                if len(key) == 0: continue
                break

            if len(key.split("_")) >= 3:
                components = key.split("_")
                trgr_nr = components[0]
                if components[0].startswith("DECOY"):
                    trgr_nr = components[1]
                try:
                    trgr_nr = int(trgr_nr)
                    return True
                except ValueError:
                    print "Format determination: Could not convert", trgr_nr, "to int."
                    return False

        return False

    def _group_precursors_by_sequence(self):
        """Group together precursors with the same charge state"""
        self._sequences_mapping = {}
        for precursor in self._precursor_mapping.keys():
            seq = precursor.split("/")[0]
            tmp = self._sequences_mapping.get(seq, [])
            tmp.append(precursor)
            self._sequences_mapping[seq] = tmp

    #
    ## Getters (data) -> see ChromatogramTransition.getData
    #
    def get_data_for_transition(self, chrom_id):
        c = self._run[str(chrom_id)] 
        return [ [c.time, c.i] ]

    def getTransitionCount(self):
        if self._in_memory:
            return len(self._run)
        else:
            return len(self._run.info['offsets']) -2

    def get_data_for_precursor(self, precursor):
        """Retrieve data for a specific precursor - data will be as list of
        pairs (timearray, intensityarray)"""

        if not self._precursor_mapping.has_key(str(precursor)):
            return [ [ [0], [0] ] ]

        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            chromatogram = self._run[str(chrom_id)] 
            transitions.append([chromatogram.time, chromatogram.i])

        if len(transitions) == 0: 
            return [ [ [0], [0] ] ]

        return transitions

    def get_range_data(self, precursor):
        return self._range_mapping.get(precursor, [0,0])

    def get_score_data(self, precursor):
        return self._score_mapping.get(precursor, None)

    def get_intensity_data(self, precursor):
        return self._intensity_mapping.get(precursor, None)

    def get_id(self):
        return self._basename

    #
    ## Getters (info)
    #
    def get_transitions_for_precursor(self, precursor):
        """
        Return the transition names for a specific precursor
        """
        return self._precursor_mapping.get(str(precursor), [])

    def get_transitions_with_mass_for_precursor(self, precursor):
        """
        Return the transition names prepended with the mass for a specific precursor
        """
        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            chromatogram = self._run[str(chrom_id)] 
            mz = chromatogram['precursors'][0]['mz']
            transitions.append(str(mz) + " m/z (" + chrom_id + ")")
        return transitions

    def get_precursors_for_sequence(self, sequence):
        return self._sequences_mapping.get(sequence, [])

    def get_all_precursor_ids(self):
        return self._precursor_mapping.keys()

    def get_all_peptide_sequences(self):
        return self._sequences_mapping.keys()

class SwathRun(object):
    """Data model for an individual SWATH injection, may contain multiple mzML files.

    This contains the model for all data from a single run (e.g. one panel in
    the viewer - in reality this could be multiple actual MS runs since in SRM
    not all peptides can be measured in the same injection or just multiple
    files generated by SWATH MS.
    """
    def __init__(self, files, runid=None):
        self.all_swathes = {}
        self.runid = runid
        self._in_memory = False
        self._files = files

        # extra info
        self._range_mapping = {}
        self._score_mapping = {}
        self._intensity_mapping = {}

        self._loadFiles(files)
        self._initialize()

    def _loadFiles(self, files):
        import pymzml
        for f in files:
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run_.original_file = f
            first = run_.next()
            mz = first['precursors'][0]['mz']
            self.all_swathes[ int(mz) ] = RunDataModel(run_, f)

    def _initialize(self):
        """ A precursor can be mapped uniquely to a certain SWATH window.
        """
        self._precursor_run_map = {}
        self._sequences_mapping = {}
        self._chrom_id_run_map = {}
        for run_key, run in self.all_swathes.iteritems():
            for key in run._precursor_mapping:
                self._precursor_run_map[key] = run_key
            for key in run._sequences_mapping:
                tmp = self._sequences_mapping.get(key, [])
                tmp.extend( run._sequences_mapping[key] )
                self._sequences_mapping[key] = tmp

    def remove_precursors(self, toremove):
        for run_key, run in self.all_swathes.iteritems():
            for key in toremove:
                run._precursor_mapping.pop(key, None)
                self._precursor_run_map.pop(key, None)
            run._group_precursors_by_sequence()

        # Re-initialize self to produce correct mapping
        self._initialize()

    #
    ## Getters (info)
    #
    def get_precursors_for_sequence(self, sequence):
        return self._sequences_mapping.get(sequence, [])

    def get_transitions_for_precursor(self, precursor):
        run = self._precursor_run_map.get( str(precursor), None)
        if run is None:
            return []
        return self.all_swathes[run].get_transitions_for_precursor(precursor)

    def get_transitions_for_precursor_display(self, precursor):
        run = self._precursor_run_map.get( str(precursor), None)
        if run is None:
            return []
        return self.all_swathes[run].get_transitions_with_mass_for_precursor(precursor)

    def get_all_precursor_ids(self):
        return self._precursor_run_map.keys()

    def get_all_peptide_sequences(self):
        res = set([])
        for m in self.all_swathes.values():
            res.update( m._sequences_mapping.keys() )
        return res

    def getAllSwathes(self):
        return self.all_swathes

    #
    ## Getters (data) -> see ChromatogramTransition.getData
    #
    def get_data_for_transition(self, chrom_id):
        raise Exception("Not implemented")

    def getTransitionCount(self):
        return sum([r.getTransitionCount() for r in self.all_swathes.values()] )

    def get_data_for_precursor(self, precursor):
        run = self._precursor_run_map[str(precursor)]
        return self.all_swathes[run].get_data_for_precursor(precursor)

    def get_range_data(self, precursor):
        return self._range_mapping.get(precursor, [0,0])

    def get_score_data(self, precursor):
        return self._score_mapping.get(precursor, None)

    def get_intensity_data(self, precursor):
        return self._intensity_mapping.get(precursor, None)

    def get_id(self):
        fileid = ""
        if len(self._files) > 0:
            fileid = os.path.basename(self._files[0]) 

        return self.runid + fileid

class SwathRunCollection(object):
    """A collection of SWATH files

    Contains multiple SwathRun objects which each represent one injection.
    """

    def __init__(self):
        self.swath_chromatograms = {}

    def initialize_from_directories(self, runid_mapping):
        """Initialize from a directory, assuming that all .mzML files in the
        same directory are from the same run.

        Requires a dictionary of form { run_id : directory }
        """
        self.swath_chromatograms = {}
        for runid, dname in runid_mapping.iteritems():
            import glob
            files = glob.glob(os.path.join(dname + "/*.mzML") )
            self.swath_chromatograms[ runid ] = SwathRun(files, runid)

    def initialize_from_chromatograms(self, runid_mapping):
        """Initialize from a set of chromatograms.

        Requires a dictionary of form { run_id : [chromatogram_files] }
        """
        self.swath_chromatograms = {}
        for runid, chromfiles in runid_mapping.iteritems():
            self.swath_chromatograms[ runid ] = SwathRun(chromfiles, runid)

    def initialize_from_files(self, filenames):
        """Initialize from individual files, setting the runid as increasing
        integers.  
        """
        self.swath_chromatograms = {}
        for i,f in enumerate(filenames):
            runid = i
            self.swath_chromatograms[ runid ] = SwathRun([f], str(runid) )

    def getSwathFiles(self):
        return self.swath_chromatograms.values()

    def getSwathFile(self, key):
        return self.swath_chromatograms[key]

    def getRunIds(self):
        return self.swath_chromatograms.keys()

class PrecursorModel():

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

    def __init__(self):
        self.runs = []

    def loadFiles(self, filenames):

        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_files(filenames)
        self.runs = [run for run in swathfiles.getSwathFiles()]

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

    def get_runs(self):
        return self.runs


    def loadFiles_with_peakgroups(self, RawData, aligned_pg_files):

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        try:
            swathfiles.initialize_from_directories( dict( [ (d["id"], d["directory"]) for d in RawData] ) )
        except KeyError:
            swathfiles.initialize_from_chromatograms( dict( [ (d["id"], d["chromatograms"]) for d in RawData] ) )
        self.runs = [run for run in swathfiles.getSwathFiles()]
        print "Find in total a collection of %s runs." % len(swathfiles.getRunIds() )

        try:
            self.read_trafo(RawData)
        except IOError:
            self.read_peakgroup_files(aligned_pg_files, swathfiles)

    def loadMixedFiles(self, rawdata_files, aligned_pg_files):
        """ Load files that contain raw data files and aligned peakgroup files.

        Since no mapping is present here, we need to infer it from the data.
        Basically, we try to map the column align_runid to the filenames of the
        input .chrom.mzML hopeing that the user did not change the filenames.
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
                raise Exception("need column header align_origfilename and align_runid")

            for this_row in reader:

                # 1. Get the original filename (find a non-NA entry) and the corresponding run id
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
        self.read_peakgroup_files(aligned_pg_files, swathfiles)
        print "Find in total a collection of %s runs." % len(swathfiles.getRunIds() )
        

    def read_trafo(self, trafo_filenames):
        # Read the transformations
        transformation_collection_ = TransformationCollection()
        for filename in [d["trafo_file"] for d in trafo_filenames]:
          transformation_collection_.readTransformationData(filename)
        transformation_collection_.initialize_from_data(reverse=True)

    def read_peakgroup_files(self, aligned_pg_files, swathfiles):
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
                    
    def load_from_yaml(self, yamlfile):

        import yaml
        data = yaml.load(open(yamlfile) )["AlignedSwathRuns"]
        alignment_files = data["PeakGroupData"]
        trafo_fnames = [d["trafo_file"] for d in data["RawData"]]
        self.loadFiles_with_peakgroups(data["RawData"], alignment_files)


CHROMTYPES = {
    0 : "Peptide", 
    1 : "Precursor", 
    2 : "Transition"
} 

CHROMTYPES_r = dict([ (v,k) for k,v in CHROMTYPES.iteritems()])

class ChromatogramTransition(object): # your internal structure

    def __init__(self, name, charge, subelements, peptideSequence=None, fullName=None, datatype="Precursor"):
        self._name = name
        self._charge = charge
        self._fullName = fullName
        self._peptideSequence = peptideSequence
        self._subelements = subelements
        self.mytype = CHROMTYPES_r[datatype]

    def getSubelements(self):
        return self._subelements

    def getPeptideSequence(self):
        if self._peptideSequence is None:
            return self.getName()
        return self._peptideSequence

    def getName(self):
        return self._name

    def getCharge(self):
        return self._charge

    def getType(self):
        return CHROMTYPES[self.mytype]

    def getData(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_data_for_precursor(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_data_for_precursor(prec[0]) 
            else:
                final_data = []
                # Sum up the data for all individual precursors
                for p in prec:
                    timedata = None
                    intdata = None
                    import numpy
                    for data in run.get_data_for_precursor(p):
                        if timedata is None:
                            timedata = numpy.array(data[0])
                            intdata = numpy.array(data[1])
                        else:
                            intdata = intdata + numpy.array(data[1])
                    final_data.append( [timedata, intdata] )
                return final_data
        elif CHROMTYPES[self.mytype] == "Transition" :
            return run.get_data_for_transition(self.getName()) 
        return [ [ [0], [0] ] ]

    def getRange(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_range_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_range_data(prec[0]) 
        return [ 0,0]

    def getProbScore(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_score_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_score_data(prec[0]) 
            else: return None
        return 1.0

    def getIntensity(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_intensity_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_intensity_data(prec[0]) 
            else: return None
        return -1.0

    def getLabel(self, run):
        """
        Get the labels for a curve (corresponding to the raw data from getData
        call) for a certain object.

        If we have a single precursors or a peptide with only one precursor, we
        show the same data as for the precursor itself. For a peptide with
        multiple precusors, we show all precursors as individual curves. For a
        single transition, we simply plot that transition.
        """
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_transitions_for_precursor_display(self.getName())
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_transitions_for_precursor_display(prec[0])
            else:
                # Peptide view with multiple precursors
                return prec
        elif CHROMTYPES[self.mytype] == "Transition" :
            return [self.getName()]
        return [ "" ]

