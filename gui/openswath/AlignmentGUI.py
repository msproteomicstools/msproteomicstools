
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

from random import random

from PyQt4 import QtGui, QtCore

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
        #print "initialize, has key", self._run.info['offsets'].has_key("DECOY_59948_YNFSDFKIPLVGNTEANIM[147]EK/3_y3")
        self.scan_for_precursor()

    def scan_for_precursor(self):
        for chrom in self._run:
            if chrom.has_key('precursors'):
                # print chrom['precursors']
                if len(chrom['precursors']) > 0:
                    if chrom['precursors'][0]['userParams'].has_key("peptide_sequence"):
                        this_prec = chrom['precursors'][0]['userParams']["peptide_sequence"] 
                        r = self._precursor_mapping.get(this_prec, [])
                        r.append(chrom['id'])
                        self._precursor_mapping[this_prec] = r

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

    def get_all_display_ids(self):
        # return self._run.info['offsets'].keys() 
        return self._precursor_mapping.keys()

class DataModel():

    def __init__(self):
        self.precursors = set([])

    def loadFiles(self, filenames):

        self.runs = []
        for f in filenames:
            print "read file", f
            import pymzml
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run = RunDataModel(run_, f)
            self.runs.append(run)
            self.precursors.update(run.get_all_display_ids())
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

    def get_runs(self):
        return self.runs

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# every model needs to have
#  - rowCount(self, parent)
#  - data(self, index, role)

class ExamplePeptides_(QtCore.QAbstractListModel):
    def rowCount(self, parent):
        return 5

    def setHorizontalHeaderLabels(self, parent):
        return 5

    def data(self, index, role):
        # DecorationRole
        # ToolTipRole
        if role == QtCore.Qt.DisplayRole:
            if index.row() == 0:
                return "test"
            elif index.row() == 1:
                return "test middle"
            else:
                return "test_later"

        if role == QtCore.Qt.ToolTipRole:
            if index.row() == 0:
                return "test Tipp"
            elif index.row() == 1:
                return "test Tipp middle"
            else:
                return "Tipp"

    def headerData(self, section, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            return "header!"

    def columnCount(self, parent):
        return 1

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

class ExamplePeptidesTreeView( QtGui.QTreeView ):

    def __init__(self):
        super(ExamplePeptidesTreeView, self).__init__()

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.openMenu)
        
    def openMenu(self, position):
    
        indexes = self.selectedIndexes()
        if len(indexes) > 0:
        
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        
        # QWidget.tr == translate function
        menu = QtGui.QMenu()
        if level == 0:
            menu.addAction(self.tr("Edit person"))
        elif level == 1:
            menu.addAction(self.tr("Edit object/container"))
        elif level == 2:
            menu.addAction(self.tr("Edit object"))
        
        menu.exec_(self.viewport().mapToGlobal(position))

class CurvePlotView(CurvePlot):

    def __init__(self, parent):
        super(CurvePlotView, self).__init__(parent)
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

    def create_curves(self, nr):

        for i in range(nr):
            param = CurveParam()
            param.label = 'My curve'
            print "try to get colors", i
            if i >= len(self.colors):
                color = COLORS.get(self.colors[0],  self.colors[0] )
            else:
                color = COLORS.get(self.colors[i],  self.colors[i] )
            param.line.color = color

            # create a new curve
            curve = CurveItemModel(param)
            self.curves.append(curve)
            self.add_item( curve )

    def update_all_curves(self, data):
        for curve in self.curves:
            curve.set_data( range( 0, 20, 2), map( lambda _: random(), range( 0, 10 ) ) )

        self.replot( )

class CurveDialogView(CurveDialog):

    def __init__(self, *args, **kwargs):
        super(CurveDialogView, self).__init__(*args, **kwargs)
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

    def update_all_curves(self, precursor):

        data = self.run.get_data_for_precursor(precursor) 
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

class CurveItemModel(CurveItem):

    def __init__(self, *args, **kwargs):
        super(CurveItemModel, self).__init__(*args, **kwargs)

        # self.initialize()

    def initialize():
        # self.curve = make.curve( [ ], [ ], "curve1", QtGui.QColor( 255, 0, 0) )
        pass

# The widget for the Graphing area on the right
class GraphArea(QtGui.QWidget):
    def __init__(self):
        super(GraphArea, self).__init__()

        self.initUI()
        self._wcount = 1
        self.c = Communicate()
        self.c.catch_mouse_press.connect(self.react_to_mouse)       
        self.c.catch_mouse_release.connect(self.react_to_mouse_release)
        
    def set_communicate(self, comm):
        self.c = comm

    def react_to_mouse(self):
        # print "react to mouse"
        pass

    def react_to_mouse_release(self):
        # print "react to release mouse"
        pass

    def initUI(self):
        self.layout = QtGui.QGridLayout(self)

    def delete_all(self):
        for i in range(self.layout.count()):
            self.layout.itemAt(i).widget().close()
        self._wcount = 1

    def add_new(self, l):
        self.layout.addWidget(l, self._wcount, 0)
        self._wcount += 1

    def mousePressEvent(self, event):
        
        self.c.catch_mouse_press.emit()

    def mouseReleaseEvent(self, event):
        
        self.c.catch_mouse_release.emit()

class ApplicationController:
    
    def __init__(self):
        pass
        
    def start(self):

        pass
    
class ApplicationView(QtGui.QWidget):
    
    def __init__(self):
        super(ApplicationView, self).__init__()
        
        self.initUI()
        
    def initUI(self):

        # self._precursor_model = ExamplePeptides_()
        self._precursor_model = ExamplePeptidesModel()
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

        # connect the two
        # self.treeView.clicked.connect(self.treeViewClicked)
        # self.treeView.selectionChanged.connect(self.treeViewClicked)
        ### self.treeView.connect(self,  QtCore.SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),  self.treeViewClicked) 
        self.treeView.selectionModel().selectionChanged.connect(self.treeViewClicked) 

        print self.treeView.selectionModel()

        self.add_plots_dummy()

    def get_precursor_model(self):
        return self._precursor_model

    def set_communication(self, c):
        self.c = c

    def add_plots_dummy(self):
        
        self.plots = []

        self.plot = CurveDialogView(edit=False, toolbar=False )
        self.plot.create_curves(3)
        self.graph_layout.add_new(self.plot)
        self.plots.append(self.plot)

        self.plot2 = CurvePlotView( self )
        self.plot2.create_curves(2)
        self.graph_layout.add_new(self.plot2)
        self.plots.append(self.plot2)

    def add_plots(self, datamodel):
        
        self.plots = []
        self.graph_layout.delete_all()

        for run in datamodel.get_runs():

            self.plot = CurveDialogView(edit=False, toolbar=False)
            self.plot.setDataModel(run)
            self.graph_layout.add_new(self.plot)
            self.plots.append(self.plot)

    def treeViewClicked(self, newvalue, oldvalue):

        assert len(self.treeView.selectedIndexes()) == 1
        selected_precursor = self.treeView.selectedIndexes()[0].data().toPyObject()

        for pl in self.plots:
            pl.update_all_curves(selected_precursor)

    def widgetclicked(self, value):
        print "clicked iittt"


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

        print "testst"
        
        # Load the files
        self.data_model.loadFiles(pyFileList)

        # get precursors from data and set it 
        pr_list = self.data_model.get_precursor_list()
        self.application.get_precursor_model().set_precursor_data(pr_list)

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
