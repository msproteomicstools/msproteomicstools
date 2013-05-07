#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from AlignmentGUI import *
from models.MSData import RunDataModel

from msproteomicstoolslib.format.TransformationCollection import TransformationCollection
from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader
from msproteomicstoolslib.algorithms.alignment.AlignmentHelper import AlignmentExperiment as Experiment 

# options
realign_runs = True
fdr_cutoff = 0.01
ONLY_SHOW_QUANTIFIED = False
ONLY_SHOW_QUANTIFIED = True

class SwathRun(object):
    def __init__(self, files):
        self.all_swathes = {}
        self.files = files
        for f in self.files:
            import pymzml
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run_.original_file = f
            first = run_.next()
            mz = first['precursors'][0]['mz']
            self.all_swathes[ int(mz) ] = run_

    def getAllSwathes(self):
        return self.all_swathes

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
            self.swath_chromatograms[ runid ] = SwathRun(files)

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


class DataModelNew(DataModel):

    def __init__(self):
        super(DataModelNew, self).__init__()
        pass

    def loadFiles(self, trafo_filenames, aligned_pg_files):

        # load new files, clean up ...
        self.runs = []
        self.precursors = set([])

        reader = SWATHScoringReader.newReader(aligned_pg_files, "openswath", readmethod="complete")
        new_exp = Experiment()
        new_exp.runs = reader.parse_files(realign_runs)
        multipeptides = new_exp.get_all_multipeptides(fdr_cutoff, verbose=False)

        # Read the transformations
        transformation_collection_ = TransformationCollection()
        for filename in [d["trafo_file"] for d in trafo_filenames]:
          transformation_collection_.readTransformationData(filename)
        transformation_collection_.initialize_from_data(reverse=True)

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_directories( dict( [ (d["id"], d["directory"]) for d in trafo_filenames] ) )
        for runid,fileobj in swathfiles.getSwathFiles().iteritems():
            for mz,run_ in fileobj.getAllSwathes().iteritems():
                run = RunDataModel(run_, run_.original_file)
                self.precursors.update(run.get_all_precursor_ids())
                run.runid = runid
                self.runs.append(run)

        # map the read-in peakgroups
        peakgroup_map = {}
        for m in multipeptides:
            pg = m.find_best_peptide_pg()
            peakgroup_map[ pg.get_value("FullPeptideName") + "/" + pg.get_value("Charge")] = m

        for rundatamodel in self.runs:
            print rundatamodel
            if ONLY_SHOW_QUANTIFIED:
                intersection = set(rundatamodel._precursor_mapping.keys()).intersection( peakgroup_map.keys() )
                rundatamodel._precursor_mapping = dict( [(k,rundatamodel._precursor_mapping[k]) for k in intersection] )
            rundatamodel._group_precursors_by_sequence()
            for key in rundatamodel._precursor_mapping.keys():
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
        self.loadFiles(data["RawData"], alignment_files)

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



