
#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ZetCode PyQt4 tutorial 

In this example, we connect a signal
of a QtGui.QSlider to a slot 
of a QtGui.QLCDNumber. 

author: Jan Bodnar
website: zetcode.com 
last edited: October 2011
"""

import sys

import re
from random import random

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, QModelIndex

from PyQt4 import Qwt5

from guiqwt.plot import CurvePlot
from guiqwt.curve import CurveItem
from guiqwt.builder import make
from guiqwt.styles import CurveParam

from guiqwt.plot import CurveDialog

from guiqwt.styles import COLORS


class Communicate(QtCore.QObject):
    
    catch_mouse_press = QtCore.pyqtSignal() 
    catch_mouse_release = QtCore.pyqtSignal() 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class RunDataModel():

    def __init__(self, run, filename):
        self._run = run
        self._filename = filename
        self._precursor_mapping = {}
        self._sequences_mapping = {}
        #print "initialize, has key", self._run.info['offsets'].has_key("DECOY_59948_YNFSDFKIPLVGNTEANIM[147]EK/3_y3")

        self._scan_for_precursor()
        self._build_sequence_mapping()

    def _scan_for_precursor(self):
        """
        Populate the mapping between precursors and the chromatogram ids.

        The precursor is of type 'PEPT[xx]IDE/3' 
        """
        openswath_format = False
        if len( self._run.info['offsets'] ) > 0:
            keys = self._run.info['offsets'].keys()
            for key in self._run.info['offsets'].keys():
                if key in ("indexList", "TIC"): continue
                break

            if len(key.split("_")) in [3,4]:
                components = key.split("_")
                trgr_nr = components[0]
                if components[0].startswith("DECOY"):
                    trgr_nr = components[1]
                try:
                    trgr_nr = int(trgr_nr)
                    openswath_format = True
                except ValueError:
                    openswath_format = False

        if openswath_format:
            if len( self._run.info['offsets'] ) > 0:
                for key in self._run.info['offsets'].keys():

                    components = key.split("_")
                    if key in ("indexList", "TIC"): continue

                    trgr_nr = str(components[1])
                    if components[0].startswith("DECOY"):
                        trgr_nr = str(components[2])

                    if self._precursor_mapping.has_key(trgr_nr):
                        self._precursor_mapping[trgr_nr].append(key)
                    else:
                        self._precursor_mapping[trgr_nr] = [key]
        else:
            # TODO fallback option!!!
            pass

    def _scan_for_precursor_by_peptide_seq(self):
        # TODO group by id if it is present and in the correct format!!!
        for chrom in self._run:
            if chrom.has_key('precursors'):
                # print chrom['precursors']
                if len(chrom['precursors']) > 0:
                    if chrom['precursors'][0]['userParams'].has_key("peptide_sequence"):
                        this_prec = chrom['precursors'][0]['userParams']["peptide_sequence"] 
                        r = self._precursor_mapping.get(this_prec, [])
                        r.append(chrom['id'])
                        self._precursor_mapping[this_prec] = r

    def get_transitions_for_precursor(self, precursor):
        return self._precursor_mapping.get(str(precursor), [])

    def get_data_for_transition(self, chrom_id):
        c = self._run[str(chrom_id)] 
        return [ [c.time, c.i] ]

    def get_data_for_precursor(self, precursor):

        # print "will try to get", precursor, "from", self._run, "at", self._filename
        # print "initialize, has key", self._run.info['offsets'].has_key("DECOY_59948_YNFSDFKIPLVGNTEANIM[147]EK/3_y3")
        # print "is equal", "DECOY_59948_YNFSDFKIPLVGNTEANIM[147]EK/3_y3" == precursor
        # print "is equal", "DECOY_59948_YNFSDFKIPLVGNTEANIM[147]EK/3_y3" == str(precursor)
        # print "prec mapping", self._precursor_mapping

        if not self._precursor_mapping.has_key(str(precursor)):
            # print "no precursor mapping"
            return [ [ [0], [0] ] ]

        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            # print "will try to get with", chrom_id
            c = self._run[str(chrom_id)] 
            # print c['id']
            transitions.append([c.time, c.i])

        #current_chroms = self._run[str(precursor)]
        # current_chroms = self._run["DECOY_59948_YNFSDFKIPLVGNTEANIM[147]EK/3_y3"]
        if len(transitions) == 0: 
            # print "transitoin len zero"
            return [ [ [0], [0] ] ]
        return transitions

    def get_precursors_for_sequence(self, sequence):
        return self._sequences_mapping.get(sequence, [])

    def get_all_precursor_ids(self):
        # return self._run.info['offsets'].keys() 
        return self._precursor_mapping.keys()

    def _build_sequence_mapping(self):
        self._sequences_mapping = {}
        for precursor in self._precursor_mapping.keys():
            seq = re.sub("[^A-Z]", "", precursor)
            tmp = self._sequences_mapping.get(seq, [])
            tmp.append(precursor)
            self._sequences_mapping[seq] = tmp
            # print "append to seq mapping", seq, tmp
        # print "seq mapping is ", self._sequences_mapping

    def get_all_peptide_sequences(self):
        return self._sequences_mapping.keys()

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

class DataModel():

    def __init__(self):
        self.precursors = set([])
        self.runs = []

    def loadFiles(self, filenames):

        self.runs = []
        for f in filenames:
            print "read file", f
            import pymzml
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run = RunDataModel(run_, f)
            self.runs.append(run)
            self.precursors.update(run.get_all_precursor_ids())
            ## first = run.next()
            ## first['product']
            ## first['precursors']
            ## mz = first['precursors'][0]['mz']
            ## # print mz
            ## all_swathes[ int(mz) ] = run
        # swath_chromatograms[ runid ] = all_swathes

        # # get the chromatograms
        # r = select_correct_swath(swath_chromatograms, current_mz)
        # chrom_ids = pg.get_value("aggr_Fragment_Annotation").split(";")
        # # print r, current_mz
        # if not r.has_key(current_run.get_id()):
        #     return ["NA"]
        # allchroms = r[current_run.get_id()]
        # current_chroms = [allchroms[chr_id] for chr_id in chrom_ids]
        # alls = 0
        # for c in current_chroms:
        #     alls += sum( [p[1] for p in c.peaks if p[0] > this_run_lwidth and p[0] < this_run_rwidth ])
        # return [alls]

    def get_precursor_list(self):
        return self.precursors

    def get_precursor_tree(self):
        # self._build_tree()
        # return self.precursors
        return self._build_tree()

    def _build_tree(self):
        # assume we have fully loaded
        peptide_sequences = set([])
        for r in self.runs:
            peptide_sequences.update( r.get_all_peptide_sequences() )

        ## The peptide sequences are our top-level items
        # print "pepseqs", peptide_sequences
        elements = []
        for seq in peptide_sequences:
            # get all precursors from all runs
            precursors = set([])
            for r in self.runs:
                precursors.update( r.get_precursors_for_sequence(seq) )
            # print "found precursros", precursors
            pelements = []
            for p in precursors:
                # get all transitions from all runs
                transitions = set([])
                for r in self.runs:
                    transitions.update( r.get_transitions_for_precursor(p) )
                tr_elements = []
                pm = PrecursorModel(p)
                for tr in transitions:
                    tr_elements.append(ChromatogramTransition(tr, -1, [], fullName=tr,
                       peptideSequence = pm.getFullSequence(), datatype="Transition") )
                pelements.append(ChromatogramTransition(p, pm.getCharge(), tr_elements, 
                       peptideSequence = pm.getFullSequence(), datatype="Precursor") )
            elements.append(ChromatogramTransition(seq, "NA", pelements, datatype="Peptide", 
                       peptideSequence=pm.getFullSequence()) )
        return elements

    def get_runs(self):
        return self.runs

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

## Generic tree models
# from http://www.hardcoded.net/articles/using_qtreeview_with_qabstractitemmodel.htm
from TreeModels import TreeNode
from TreeModels import TreeModel

CHROMTYPES = {
    0 : "Peptide", 
    1 : "Precursor", 
    2 : "Transition"
} 

CHROMTYPES_r = dict([ (v,k) for k,v in CHROMTYPES.iteritems()])

class ChromatogramTransition(object): # your internal structure
    def __init__(self, name, charge, subelements, peptideSequence=None, fullName=None, datatype="Precursor"):
        self.name = name
        self.charge = charge
        self.fullName = fullName
        self.peptideSequence = peptideSequence
        self.subelements = subelements
        self.mytype = CHROMTYPES_r[datatype]

    def getPeptideSequence(self):
        if self.peptideSequence is None:
            return self.name
        return self.peptideSequence

    def getName(self):
        return self.name

    def getType(self):
        return CHROMTYPES[self.mytype]

    def getData(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            # Precursor type
            return run.get_data_for_precursor(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_data_for_precursor(prec[0]) 
            else:
                # TODO dont just show the first one!
                pass
                # return run.get_data_for_precursor(prec[0]) 
        elif CHROMTYPES[self.mytype] == "Transition" :
            return run.get_data_for_transition(self.getName()) 
        return [ [ [0], [0] ] ]

# A TreeNode element
class PeptideTreeNode(TreeNode):
    def __init__(self, ref, parent, row):
        self.ref = ref
        TreeNode.__init__(self, parent, row)

    def _getChildren(self):
        return [PeptideTreeNode(elem, self, index)
            for index, elem in enumerate(self.ref.subelements)]

class PeptideTree(TreeModel):
    def __init__(self, rootElements):
        self.rootElements = rootElements
        TreeModel.__init__(self)

    def _getRootNodes(self):
        return [PeptideTreeNode(elem, None, index)
            for index, elem in enumerate(self.rootElements)]

    def columnCount(self, parent):
        return 2

    def data(self, index, role):
        if not index.isValid():
            return None
        node = index.internalPointer()
        if role == Qt.DisplayRole and index.column() == 0:
            return QtCore.QVariant( node.ref.getPeptideSequence() )
        if role == Qt.DisplayRole and index.column() == 1:
            return QtCore.QVariant( node.ref.charge )
        return None

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 0:
            return 'Peptide Sequence'
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 1:
            return 'Charge'
        return None

    def set_precursor_tree_structure(self, data):
        self.rootElements = []
        # print "Set with dataA", data
        self.rootElements = data

        # initialize super method again
        self.initialize()

    def set_precursor_data(self, data):
        self.rootElements = []
        for data_item in data:
            try:
                charge = data_item.split("/")[1].split("_")[0]
            except Exception:
                charge = 0
            self.rootElements.append(ChromatogramTransition(data_item, charge, [] ) )

        # initialize super method again
        self.initialize()

    def setHorizontalHeaderLabels(self, parent):
        pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# every model needs to have
#  - rowCount(self, parent)
#  - data(self, index, role)
#  - setData (
#  - headerData
#  - flags  # should be editable/selectable

# 
## Tree Model - there are two options
#
# class ExamplePeptides_(QtCore.QAbstractItemModel):
class ExamplePeptides_(QtCore.QAbstractListModel):

    def __init__(self):
        super(ExamplePeptides_, self).__init__()
        self.precursor_data = [ "a", "b", "c", "d", "c"]

    def rowCount(self, parent):
        return len(self.precursor_data)

    def index(self, row, column, parent):
        # what do to
        pass

    def setHorizontalHeaderLabels(self, parent):
        pass

    def data(self, index, role):
        # DecorationRole
        # ToolTipRole
        if index.isValid() and role == QtCore.Qt.DisplayRole:
            #return QtCore.QVariant( self.precursor_data[ index.row() ] )
            return self.precursor_data[ index.row() ]
            

        # if role == QtCore.Qt.ToolTipRole:
        #     if index.row() == 0:
        #         return "test Tipp"
        #     elif index.row() == 1:
        #         return "test Tipp middle"
        #     else:
        #         return "Tipp"

    def headerData(self, section, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            return "header!"

    def columnCount(self, parent):
        return 1

    def set_precursor_data(self, data):
        self.precursor_data = []
        for data_item in data:
            self.precursor_data.append(data_item)

class ExamplePeptidesModel(QtGui.QStandardItemModel):

    def __init__(self):
        super(ExamplePeptidesModel, self).__init__()

        self.initialize()

    def initialize(self):
        item = QtGui.QStandardItem("root")
        self.appendRow(item)
        item.appendRow( QtGui.QStandardItem("a") )
        item.appendRow( QtGui.QStandardItem("foo") )
        item2 = QtGui.QStandardItem("root2")
        item2.appendRow( QtGui.QStandardItem("foo") )
        self.appendRow(item2)

    def set_precursor_data(self, data):
        self.clear()
        for data_item in data:
            item = QtGui.QStandardItem(data_item)
            self.appendRow(item)

# 
## Tree View
#
class ExamplePeptidesTreeView( QtGui.QTreeView ):

    def __init__(self):
        super(ExamplePeptidesTreeView, self).__init__()

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.openMenu)
        
    def openMenu(self, position):
        pass
    #
    #    indexes = self.selectedIndexes()
    #    if len(indexes) > 0:
    #    
    #        level = 0
    #        index = indexes[0]
    #        while index.parent().isValid():
    #            index = index.parent()
    #            level += 1
    #    
    #    # QWidget.tr == translate function
    #    menu = QtGui.QMenu()
    #    if level == 0:
    #        menu.addAction(self.tr("Edit person"))
    #    elif level == 1:
    #        menu.addAction(self.tr("Edit object/container"))
    #    elif level == 2:
    #        menu.addAction(self.tr("Edit object"))
    #    
    #    menu.exec_(self.viewport().mapToGlobal(position))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# A single curve
class CurveItemModel(CurveItem):

    def __init__(self, *args, **kwargs):
        super(CurveItemModel, self).__init__(*args, **kwargs)

# 
## The widget for a single plot on the right
#
class MultiLinePlot(CurveDialog):
    """For the Curve window we could use a CurveDialog or a CurvePlot. 

    CurveDialog has more features and seems more advanced.
    """

    def __init__(self, *args, **kwargs):
        super(MultiLinePlot, self).__init__(*args, **kwargs)
        self.myrange = None
        self.run = None
        self.initialize()

    def initialize(self):
        self.colors = [
                QtGui.QColor( 255, 0, 0),
                QtGui.QColor( 50, 50, 50),
                QtGui.QColor( 255, 0, 255),
                QtGui.QColor( 0, 200, 100),
                QtGui.QColor( 0, 0, 255),
                QtGui.QColor( 255, 0, 80),
                QtGui.QColor( 100, 0, 80),
                QtGui.QColor( 100, 0, 0)
        ]
        self.curves = []
        self.ranges = []

    def setDataModel(self, run):
        self.run = run

    def create_curves(self, nr):

        self.curves = []
        plot = self.get_plot()
        plot.del_all_items(except_grid=False)
        for i in range(nr):
            param = CurveParam()
            param.label = 'My curve'
            #color = COLORS.get(self.colors[i],  self.colors[i] )
            if i >= len(self.colors):
                color = COLORS.get(self.colors[0],  self.colors[0] )
            else:
                color = COLORS.get(self.colors[i],  self.colors[i] )
 
            param.line.color = color

            # create a new curve
            curve = CurveItemModel(param)

            self.curves.append(curve)

            plot.add_item( curve )

        self.myrange = make.range(0,0)
        self.myrange.itemChanged()
        # disp2 = make.computations(self.myrange, "TL",
        #                               [(curve, "min=%.5f", lambda x,y: y.min()),
        #                                (curve, "max=%.5f", lambda x,y: y.max()),
        #                                (curve, "avg=%.5f", lambda x,y: y.mean())])
        # plot.add_item( disp2 )
        plot.add_item( self.myrange )

    def rangeChanged(self, data):
        print "range changed"

    def update_all_curves(self, chr_transition):

        # if chr_transition.getType() == "Precursor":
        #     precursor = chr_transition.getName()

        data = chr_transition.getData(self.run) # self.run.get_data_for_precursor(precursor) 
        self.create_curves(len(data))

        for d, curve in zip(data, self.curves):
            curve.set_data( d[0], d[1] )
        ## # TODO only here we can get the range ???
        ## # print "range was ", self.myrange.get_range()
        ## r = int( random() * 15)
        ## if not self.myrange is None: self.myrange.set_range(r, r+5)

        self.get_plot().replot( )

    def mouseReleaseEvent(self, event):
        pass
        # TODO here i can capture the mouse release and the range !
        # print "mouse was released, range was ", self.myrange.get_range()

# 
## The widget for the Graphing area on the right
#
class GraphArea(QtGui.QWidget):

    def __init__(self):
        super(GraphArea, self).__init__()

        self.initUI()
        self._wcount = 1
        self.c = Communicate()
        self.plots = []
        
    def set_communicate(self, comm):
        self.c = comm

    def initUI(self):
        self.layout = QtGui.QGridLayout(self)

    def delete_all(self):
        for i in range(self.layout.count()):
            self.layout.itemAt(i).widget().close()
        self._wcount = 1

    def add_new(self, l):
        self.layout.addWidget(l, self._wcount, 0)
        self._wcount += 1

    def add_plots_dummy(self):
        
        self.plots = []

        self.plot = MultiLinePlot(edit=False, toolbar=False )
        self.plot.create_curves(3)
        self.add_new(self.plot)
        self.plots.append(self.plot)

        #self.plot2 = CurvePlotView( self )
        self.plot2 = MultiLinePlot(edit=False, toolbar=False )
        self.plot2.create_curves(2)
        self.add_new(self.plot2)
        self.plots.append(self.plot2)

    def add_plots(self, datamodel):
        
        self.plots = []
        self.delete_all()

        for run in datamodel.get_runs():

            self.plot = MultiLinePlot(edit=False, toolbar=False)
            self.plot.setDataModel(run)
            self.add_new(self.plot)
            self.plots.append(self.plot)

    def update_all_plots(self, chr_transition):
    
        for pl in self.plots:
            pl.update_all_curves(chr_transition)

class ApplicationView(QtGui.QWidget):
    
    def __init__(self):
        super(ApplicationView, self).__init__()
        
        self.initUI()
        
    def initUI(self):

        # self._precursor_model = ExamplePeptides_()
        # self._precursor_model = ExamplePeptidesModel()
        self._precursor_model = PeptideTree([])
        self._precursor_model.setHorizontalHeaderLabels([self.tr("Peptides")])

        self.treeView = ExamplePeptidesTreeView()
        self.treeView.setModel(self._precursor_model)

        self.graph_layout = GraphArea()
        horizontal_splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        horizontal_splitter.addWidget(self.treeView)
        horizontal_splitter.addWidget(self.graph_layout)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(horizontal_splitter)

        self.setLayout(hbox)

        # QItemSelectionModel -> conect the tree to here
        self.treeView.selectionModel().selectionChanged.connect(self.treeViewClicked) 

        # add dummy plots to the graph layout
        self.graph_layout.add_plots_dummy()

    def get_precursor_model(self):
        return self._precursor_model

    def set_communication(self, c):
        self.c = c

    def treeViewClicked(self, newvalue, oldvalue):

        # assert that only one single element was selected (even if multiple
        # columns are present) <=> more than one needs to be selected and they
        # all need to have the same internal pointer object
        assert len(newvalue.indexes()) > 0 
        assert all(x.internalPointer() == newvalue.indexes()[0].internalPointer() for x in newvalue.indexes())

        # selected_precursor = newvalue.indexes()[0].internalPointer().ref.getName()

        self.graph_layout.update_all_plots(newvalue.indexes()[0].internalPointer().ref)

    def add_plots(self, datamodel):
        self.graph_layout.add_plots(datamodel)

    def widgetclicked(self, value):
        print "clicked iittt"


#
## Main Window
# 
class MainWindow(QtGui.QMainWindow):
    
    def __init__(self):
        super(MainWindow, self).__init__()
        
        self.c = Communicate()
        self.data_model = DataModel()

        self.initUI()
        
    def initUI(self):               
        
        self.application = ApplicationView()
        self.application.set_communication(self.c)
        self.setCentralWidget(self.application)

        ###################################
        # Actions
        ###################################
        icon = QtGui.QIcon("web.png")

        openFile = QtGui.QAction(icon, 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.showDialog)

        exitAction = QtGui.QAction(icon, 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

        self.statusBar()

        ###################################
        # Menubar
        ###################################
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)       
        fileMenu.addAction(exitAction)

        ###################################
        # Toolbar
        ###################################
        toolbar = self.addToolBar('Exit')
        toolbar.addAction(openFile)
        # toolbar.addAction(exitAction)
        
        # self.setGeometry(300, 300, 250, 150)
        self.resize(850, 550)
        self.center()
        self.setWindowTitle('Hannes example')
        self.show()
        self.statusBar().showMessage('Ready')

        # for testing only TODO
        self.showDialog()

    def showDialog(self):

        ## fileList = QtGui.QFileDialog.getOpenFileNames(self, 'Open file', 
        ##                                               "/home/hr/projects/msproteomicstools/mzmls" )
        ## 
        ## print "opened file ", fileList
        ## pyFileList = [str(f) for f in fileList]

        pyFileList = ['/home/hr/projects/msproteomicstools/mzmls/split_hroest_K120808_Strep10PlasmaBiolRepl1_R02_SW-Strep_10%_Plasma_Biol_Repl1_16._chrom.mzML',
        '/home/hr/projects/msproteomicstools/mzmls/split_hroest_K120808_Strep10PlasmaBiolRepl1_R01_SW-Strep_10%_Plasma_Biol_Repl1_16._chrom.mzML']
        pyFileList = ['/tmp/test.mzML']

        print "testst"
        
        # Load the files
        self.data_model.loadFiles(pyFileList)

        # get precursors from data and set it 
        pr_list = self.data_model.get_precursor_list()
        precursor_model = self.application.get_precursor_model()
        if "set_precursor_tree_structure" in dir(precursor_model):
            precursor_model.set_precursor_tree_structure(self.data_model.get_precursor_tree())
        else:
            precursor_model.set_precursor_data(pr_list)

        self.application.add_plots(self.data_model)

        # print "will open files", [str(f) for f in fileList]


        # TODO open those files and process them

        # f = open(fname, 'r')
        # 
        # with f:        
        #     data = f.read()
        #     self.textEdit.setText(data) 

    def center(self):
        
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())
