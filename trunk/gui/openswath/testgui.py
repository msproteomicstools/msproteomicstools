#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from AlignmentGUI import *
from models.MSData import RunDataModel

# options
testoutfile = "/home/hr/projects/msproteomicstools/testout_all.csv"
testoutfile = "/tmp/testout.csv"
realign_runs = True
fdr_cutoff = 0.01


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

    def initialize_from_trafofiles(self, trafo_filenames):
        res  = {}
        for filename in trafo_filenames:
            # get the run id
            f = open(filename, "r")
            header = f.next().split("\t")
            f.close()
            runid = header[1]
            dname = os.path.dirname(filename)
            res[runid] = dname
        self.initialize_from_directories(res)

    def initialize_from_directories(self, runid_mapping):
        self.swath_chromatograms = {}
        for runid, dname in runid_mapping.iteritems():
            import glob
            files = glob.glob(os.path.join(dname + "/*.mzML") )
            self.swath_chromatograms[ runid ] = SwathRun(files)

    def initialize_from_files(self, filenames):
        self.swath_chromatograms = {}
        for i,f in enumerate(filenames):
            self.swath_chromatograms[ runid ] = SwathRun([f])

    def getSwathFiles(self):
        return self.swath_chromatograms


class DataModelNew(DataModel):

    def __init__(self):
        super(DataModelNew, self).__init__()
        pass

    def loadFiles(self, trafo_filenames):

        # load new files, clean up ...
        self.runs = []
        self.precursors = set([])

        sys.path.append("/home/hr/projects/msproteomicstools")
        from feature_alignment import Experiment, TransformationCollection
        from msproteomicstoolslib.format.SWATHScoringReader import SWATHScoringReader

        reader = SWATHScoringReader.newReader([testoutfile], "openswath", readmethod="complete")
        new_exp = Experiment()
        new_exp.runs = reader.parse_files(realign_runs)
        multipeptides = new_exp.get_all_multipeptides(fdr_cutoff, verbose=False)

        # Read the transformations
        transformation_collection_ = TransformationCollection()
        for filename in trafo_filenames:
          transformation_collection_.readTransformationData(filename)
        transformation_collection_.initialize_from_data(reverse=True)

        # Read the chromatograms
        swathfiles = SwathRunCollection()
        swathfiles.initialize_from_trafofiles(trafo_filenames)
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
            intersection = set(rundatamodel._precursor_mapping.keys()).intersection( peakgroup_map.keys() )
            rundatamodel._precursor_mapping = dict( [(k,rundatamodel._precursor_mapping[k]) for k in intersection] )
            rundatamodel._group_precursors_by_sequence()
            for key in rundatamodel._precursor_mapping.keys():
                m = peakgroup_map[ key ]
                if m.has_peptide(rundatamodel.runid):
                    pg = m.get_peptide(rundatamodel.runid).get_best_peakgroup()
                    try:
                        rundatamodel._range_mapping[key] = [ float(pg.get_value("leftWidth")), float(pg.get_value("rightWidth")) ]
                        rundatamodel._score_mapping[key] = float(pg.get_value("m_score"))
                        rundatamodel._intensity_mapping[key] = float(pg.get_value("Intensity"))
                    except Exception: 
                        pass


class MainWindowNew(MainWindow):
    
    def __init__(self):
        super(MainWindowNew, self).__init__()
        
        self.c = Communicate()
        self.data_model = DataModelNew()

        self.initUI()

trafo_fnames = ['analysis/alignment/strep_align/Strep0_Repl1_R02/transformation-0_0-0_1.tr', 'analysis/alignment/strep_align/Strep0_Repl2_R02/transformation-0_1-0_1.tr', 'analysis/alignment/strep_align/Strep10_Repl1_R02/transformation-0_2-0_1.tr', 'analysis/alignment/strep_align/Strep10_Repl2_R02/transformation-0_3-0_1.tr']
trafo_fnames = ['/home/hr/projects/msproteomicstools/' + f for f in trafo_fnames]



app = QtGui.QApplication(sys.argv)
ex = MainWindowNew()
ex.data_model.loadFiles(trafo_fnames)
ex._refresh_view()
sys.exit(app.exec_())



