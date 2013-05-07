#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from AlignmentGUI import *
from models.MSData import RunDataModel

from msproteomicstoolslib.format.TransformationCollection import TransformationCollection
from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import AlignmentExperiment as Experiment 

# options
REALIGN_RUNS = True
FDR_CUTOFF = 0.01
ONLY_SHOW_QUANTIFIED = False
ONLY_SHOW_QUANTIFIED = True

class SwathRun(object):
    """Data model for an individual SWATH injection, may contain multiple mzML files.
    """
    def __init__(self, files, runid=None):
        self.all_swathes = {}
        self._files = files
        self.runid = runid
        self._in_memory = False
        for f in self._files:
            import pymzml
            print "Load file", f
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run_.original_file = f
            first = run_.next()
            mz = first['precursors'][0]['mz']
            self.all_swathes[ int(mz) ] = RunDataModel(run_, f)

        self._initialize()

        # extra info
        self._range_mapping = {}
        self._score_mapping = {}
        self._intensity_mapping = {}

    def _initialize(self):
        self._precursor_run_map = {}
        self._sequences_mapping = {}
        self._chrom_id_run_map = {}
        for run_key, run in self.all_swathes.iteritems():
            for key in run._precursor_mapping:
                self._precursor_run_map[key] = run_key
            for key in run._sequences_mapping:
                # self._sequence_run_map[key] = run_key
                tmp = self._sequences_mapping.get(key, [])
                tmp.extend( run._sequences_mapping[key] )
                self._sequences_mapping[key] = tmp

    def remove_precursors(self, toremove):
        for run_key, run in self.all_swathes.iteritems():
            for key in toremove:
                run._precursor_mapping.pop(key, None)
                self._precursor_run_map.pop(key, None)
            run._group_precursors_by_sequence()

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
        return self.runid

class SwathRunCollection(object):
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

    def initialize_from_files(self, filenames):
        """Initialize from individual files, setting the runid as increasing
        integers.  
        """
        self.swath_chromatograms = {}
        for i,f in enumerate(filenames):
            runid = i
            self.swath_chromatograms[ runid ] = SwathRun([f])

    def getSwathFiles(self):
        return self.swath_chromatograms

    def getRunIds(self):
        return self.swath_chromatograms.keys()

class DataModelNew(DataModel):

    def __init__(self):
        super(DataModelNew, self).__init__()
        pass

    def loadFiles_with_peakgroups(self, trafo_filenames, aligned_pg_files):

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_directories( dict( [ (d["id"], d["directory"]) for d in trafo_filenames] ) )
        self.runs = [run for run in swathfiles.getSwathFiles().values()]
        print "Find in total a collection of %s runs." % len(swathfiles.getRunIds() )

        self.read_trafo(trafo_filenames)
        self.map_this(aligned_pg_files, swathfiles)

    def read_trafo(self, trafo_filenames):
        # Read the transformations
        transformation_collection_ = TransformationCollection()
        for filename in [d["trafo_file"] for d in trafo_filenames]:
          transformation_collection_.readTransformationData(filename)
        transformation_collection_.initialize_from_data(reverse=True)

    def map_this(self, aligned_pg_files, swathfiles):
        # map the read-in peakgroups
        reader = SWATHScoringReader.newReader(aligned_pg_files, "openswath", readmethod="complete")
        new_exp = Experiment()
        new_exp.runs = reader.parse_files(REALIGN_RUNS)
        multipeptides = new_exp.get_all_multipeptides(FDR_CUTOFF, verbose=False)

        peakgroup_map = {}
        for m in multipeptides:
            pg = m.find_best_peptide_pg()
            peakgroup_map[ pg.get_value("FullPeptideName") + "/" + pg.get_value("Charge")] = m

        for rundatamodel in swathfiles.getSwathFiles().values():
            if ONLY_SHOW_QUANTIFIED:
                intersection = set(rundatamodel.get_all_precursor_ids()).intersection( peakgroup_map.keys() )
                todelete = set(rundatamodel.get_all_precursor_ids()).difference(intersection)
                rundatamodel.remove_precursors(todelete)
            for key in rundatamodel.get_all_precursor_ids():
                if not peakgroup_map.has_key(key): continue
                m = peakgroup_map[ key ]
                if m.has_peptide(rundatamodel.runid):
                    pg = m.get_peptide(rundatamodel.runid).get_best_peakgroup()
                    try:
                        rundatamodel._range_mapping[key] = [ float(pg.get_value("leftWidth")), float(pg.get_value("rightWidth")) ]
                        rundatamodel._score_mapping[key] = float(pg.get_value("m_score"))
                        rundatamodel._intensity_mapping[key] = float(pg.get_value("Intensity"))
                    except Exception: 
                        pass
                    
    def load_from_yaml(self, yamlfile):

        import yaml
        data = yaml.load(open(yamlfile) )["AlignedSwathRuns"]
        alignment_files = data["PeakGroupData"]
        trafo_fnames = [d["trafo_file"] for d in data["RawData"]]
        self.loadFiles_with_peakgroups(data["RawData"], alignment_files)

class MainWindowNew(MainWindow):
    
    def __init__(self):
        super(MainWindowNew, self).__init__()
        
        self.c = Communicate()
        self.data_model = DataModelNew()

        self.initUI()

"""
Sample Yaml file:

AlignedSwathRuns:
  PeakGroupData: [/path/peakgroups.csv]
  RawData:
  - {directory: /path/Strep0_Repl1_R02, id: '0_0', trafo_file:  /path/Strep0_Repl1_R02/transformation-0_0-0_1.tr}
  - {directory: /path/Strep0_Repl2_R02, id: '0_1', trafo_file:  /path/Strep0_Repl2_R02/transformation-0_1-0_1.tr}
  - {directory: /path/Strep10_Repl1_R02, id: '0_2', trafo_file: /path/Strep10_Repl1_R02/transformation-0_2-0_1.tr}
  - {directory: /path/Strep10_Repl2_R02, id: '0_3', trafo_file: /path/Strep10_Repl2_R02/transformation-0_3-0_1.tr}
  ReferenceRun: '0_1'
"""

app = QtGui.QApplication(sys.argv)
ex = MainWindowNew()
ex.data_model.load_from_yaml("/home/hr/projects/msproteomicstools/test100align.yaml")
# ex.data_model.load_from_yaml("/home/hr/projects/msproteomicstools/allAlign.yaml")
ex._refresh_view()
sys.exit(app.exec_())



