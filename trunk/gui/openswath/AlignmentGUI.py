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


"""
OpenSwath Viewer

Usage:
    middle-mouse click : auto-zoom
    middle-mouse move  : pan
    right-mouse move   : zoom

TODO: 
    use QDockWidget

Structure
- AlignmentGUI.py (main)
./models
    - ./models/TreeModels.py (contains the generic tree models)
    - ./models/PeptideTree.py (contains the specific implementation of the left side peptide tree)
    - ./models/ChromatogramTransition.py (contains the chromatogram transition abstraction which is stored in the peptide tree view)

    - ./models/SingleChromatogramFile.py (data model for a single chrom.mzML file)
    - ./models/SwathRun.py (data model for a single SWATH-MS run, may contain multiple chrom.mzML files)
    - ./models/SwathRunCollection.py (data model for a set of SWATH-MS runs)
    - ./models/MSData.py (contains the mass spectrometric data models)

./views
    - ./views/PeptideTree.py (contains the tree view implementation, derived from QtGui.QTreeView)
    - ./views/Plot.py (contains the plot view, derived from Qwt.QwtPlot for the Qwt implementation or from the GuiQwt library)

"""

import sys,time, re

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, QModelIndex

# global parameters
TITLE_FONT_SIZE = 10
AXIS_FONT_SIZE = 8
AUTOSCALE_Y_AXIS = True
# Turning off guiqwt (USE_GUIQWT=False) should be safer, it then uses the
# fallback to plain Qwt.
USE_GUIQWT = True

class Communicate(QtCore.QObject):
    
    catch_mouse_press = QtCore.pyqtSignal() 
    catch_mouse_release = QtCore.pyqtSignal() 

# 
## Data Model and Tree Model
#
from models.MSData import DataModel
from models.PeptideTree import PeptideTree

# 
## Views for the plots and the peptide tree (on the right)
#
from views.PeptideTree import PeptidesTreeView
if USE_GUIQWT:
    from views.Plot import GuiQwtMultiLinePlot as MultiLinePlot
else:
    from views.Plot import QwtMultiLinePlot as MultiLinePlot

# 
## The widget for the graphing area on the right
#
class GraphArea(QtGui.QWidget):

    def __init__(self):
        super(GraphArea, self).__init__()

        self.initUI()
        self._wcount = 1
        self.plots = []
        # self.c.catch_mouse_press.connect(self.react_to_mouse)
        # self.c.catch_mouse_release.connect(self.react_to_mouse_release)
         
    @QtCore.pyqtSlot(float, float, float, float)
    def plotZoomChanged(self, xmin, xmax, ymin, ymax):
        # after (moving the image), adjust and replot _all_ plots
        self.reset_axis_all_plots([xmin, xmax], [ymin, ymax], AUTOSCALE_Y_AXIS)

    # def react_to_mouse(self):
    #     print "react to mouse press"
    # 
    # def react_to_mouse_release(self):
    #     print "react to release mouse"
    # 
    # def mousePressEvent(self, event):
    #     self.c.catch_mouse_press.emit()
    # 
    # def mouseReleaseEvent(self, event):
    #     self.c.catch_mouse_release.emit()
        
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
        self.plot.create_curves([1,2,3], [0,0])
        self.add_new(self.plot)
        self.plots.append(self.plot)

        #self.plot2 = CurvePlotView( self )
        self.plot2 = MultiLinePlot(edit=False, toolbar=False )
        self.plot2.create_curves([1,2], [0,0])
        self.add_new(self.plot2)
        self.plots.append(self.plot2)

    def add_plots(self, datamodel):
        
        self.plots = []
        self.delete_all()

        for i, run in enumerate(datamodel.get_runs()):

            self.plot = MultiLinePlot(edit=False, toolbar=False,
                                      options=dict(xlabel="Time (s)", ylabel="Intensity") )
            self.plot.setDataModel(run)

            # set font and title of plot
            self.plot.setTitleFontSize(TITLE_FONT_SIZE)
            self.plot.setAxisFontSize(AXIS_FONT_SIZE)

            self.plot.zoomChanged.connect(self.plotZoomChanged) 

            self.layout.addWidget(self.plot, i % 3, int(i/3) )
            self.plots.append(self.plot)

    def reset_axis_all_plots(self, x_range, y_range, autoscale_y=False):
        for i, pl in enumerate(self.plots):
            pl.set_x_limits(x_range[0], x_range[1])
            if autoscale_y:
                pl.set_y_limits_auto(x_range[0], x_range[1])
            else:
                pl.set_y_limits(y_range[0], y_range[1])
            pl.replot()

    def update_all_plots(self, chr_transition):
        """
        We update the plots for all runs.

        1. Load the data from the ChromatogramTransition object
        2. Plot the data for each plot, calling update_all_curves on each plot
        """
        # Get all data, compute overall min/max
        xmins = []
        xmaxs = []
        pairs = []
        # loading the data takes about 3-4 ms per plot
        for pl in self.plots:
            data = chr_transition.getData(pl.run) 
            pairs.append( data )
            xmins.extend( [min(d[0]) for d in data] )
            xmaxs.extend( [max(d[0]) for d in data] )

        for i, pl in enumerate(self.plots):
            data = pairs[i]
            labels = chr_transition.getLabel(pl.run) 
            ranges = chr_transition.getRange(pl.run) 
            mscore = chr_transition.getProbScore(pl.run) 
            intensity = chr_transition.getIntensity(pl.run) 
            # this next command takes about 10 ms per plot with Qwt, ca 30-40 ms with GuiQwt
            pl.update_all_curves(data, labels, ranges, mscore, intensity)
            pl.set_x_limits(min(xmins),max(xmaxs))
            pl.replot()

#
## Peptide Tree Widget (left side)
# 
class PeptideTreeWidget(QtGui.QWidget):

    # Signals
    selectionChanged = QtCore.pyqtSignal(QModelIndex)

    def __init__(self):
        super(PeptideTreeWidget, self).__init__()
        self.initUI()

    def initUI(self):

        # Set up the model and the view
        self._precursor_model = PeptideTree([])
        self.treeView = PeptidesTreeView()
        self.treeView.setModel(self._precursor_model)

        ## -> this would be how to sort but it seems to be buggy
        ## self.pProxyModel = QtGui.QSortFilterProxyModel()
        ## self.pProxyModel.setSourceModel(self._precursor_model)
        ## self.treeView.setModel(self.pProxyModel)
        ## self.treeView.setSortingEnabled(True)

        searchbox_layout = QtGui.QHBoxLayout()
        self.treeLineEdit = QtGui.QLineEdit()
        self.treeComboBox = QtGui.QComboBox()
        # Populate the ComboBox
        for i in range(self._precursor_model.columnCount(self)):
            self.treeComboBox.addItem(self._precursor_model.headerData(i, Qt.Horizontal, Qt.DisplayRole))

        searchbox_layout.addWidget(self.treeLineEdit) 
        searchbox_layout.addWidget(self.treeComboBox) 

        # Combine the tree and the searchbox
        self.leftside_layout = QtGui.QVBoxLayout()
        self.leftside_layout.addWidget(self.treeView)
        self.leftside_layout.addLayout(searchbox_layout)
        self.setLayout(self.leftside_layout)

        #
        ## Connect the Signals
        #
        self.treeLineEdit.textChanged.connect(self.changedTextTest)
        self.treeLineEdit.returnPressed.connect(self.changeReturnPressedTest)
        # QItemSelectionModel -> connect the tree to here
        self.treeView.selectionModel().selectionChanged.connect(self.treeViewSelectionChanged) 

    def changeReturnPressedTest(self):

        if self.treeiter is None:
            return

        column =  self.treeComboBox.currentIndex()
        try:
            model_idx = self.treeiter.next()
        except StopIteration:
            # If we have reached the bottom, wrap around
            try:
                self.treeiter = self.generate_it(self.treeView, column, self.treeLineEdit.text())
                model_idx = self.treeiter.next()
            except StopIteration:
                # No match found
                return

        self.treeView.selectAndScrollTo(model_idx)

    # Iterator generator function to go through all indexes
    def generate_it(self, treeView, column_, text_):
        s = re.compile(str(text_), re.IGNORECASE)
        m = treeView.model()
        for model_idx in treeView.iterTopLevelElements(column_):
            display_data = m.data(model_idx, Qt.DisplayRole).toPyObject()
            if s.search(display_data):
                yield model_idx

    def changedTextTest(self, text):

        column =  self.treeComboBox.currentIndex()
        try:
            self.treeiter = self.generate_it(self.treeView, column, text)
            model_idx = self.treeiter.next()
        except StopIteration:
            # No match found
            return

        self.treeView.selectAndScrollTo(model_idx)

    @QtCore.pyqtSlot(QtGui.QItemSelectionModel, QtGui.QItemSelectionModel)
    def treeViewSelectionChanged(self, newvalue, oldvalue):
        if len(newvalue.indexes()) == 0 :
            return

        # assert that the the underlying selected element is always the same. 
        assert isinstance(newvalue, QtGui.QItemSelection)
        assert all(x.internalPointer() == newvalue.indexes()[0].internalPointer() for x in newvalue.indexes())
        self.selectionChanged.emit(newvalue.indexes()[0])

    def expandLevel(self, level):
        if level == "Peptides":
            self.treeView.expandToDepth(0)
        elif level == "Precursors":
            self.treeView.expandToDepth(1)
        elif level == "smart":
            self.treeView.expandMultiElementItems()
        else:
            self.treeView.collapseAll()

    def get_precursor_model(self):
        return self._precursor_model

#
## Main Widget
# 
class ApplicationView(QtGui.QWidget):
    
    # Signals
    plotsUpdated = QtCore.pyqtSignal(float)

    def __init__(self, parent):
        super(ApplicationView, self).__init__()
        self.parent = None
        self.treeiter = None
        self.initUI()
        
    @QtCore.pyqtSlot(QModelIndex)
    def treeSelectionChanged(self, idx):
        s = time.time()
        self.graph_layout.update_all_plots(idx.internalPointer().ref)
        self.plotsUpdated.emit(time.time()-s)

    def initUI(self):

        self.leftside = PeptideTreeWidget()
        self.leftside.selectionChanged.connect(self.treeSelectionChanged)

        # Do the main application (leftside/graphing area)
        self.graph_layout = GraphArea()
        horizontal_splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        horizontal_splitter.addWidget(self.leftside)
        horizontal_splitter.addWidget(self.graph_layout)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(horizontal_splitter)

        self.setLayout(hbox)

        # add dummy plots to the graph layout
        self.graph_layout.add_plots_dummy()

    def get_precursor_model(self):
        return self.leftside.get_precursor_model()

    def set_communication(self, c):
        self.c = c


    def add_plots(self, datamodel):
        self.graph_layout.add_plots(datamodel)
        self.leftside.expandLevel("smart")

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
        
        self.application = ApplicationView(self)
        self.application.set_communication(self.c)
        self.application.plotsUpdated.connect(self.plotsUpdated)
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
        self.setWindowTitle('OpenSWATH Alignment GUI')
        self.show()
        self.statusBar().showMessage('Ready')

    @QtCore.pyqtSlot(float)
    def plotsUpdated(self, time_taken):
        self.statusBar().showMessage(self.data_model.getStatus() + ". Drawn plots in %0.4fs."  % (time_taken))

    def showDialog(self):
        """
        Show the open file dialog to load a new file.
        """

        fileList = QtGui.QFileDialog.getOpenFileNames(self, 'Open file')
        pyFileList = [str(f) for f in fileList]

        # Load the files
        start = time.time() 
        if len(pyFileList) == 1 and pyFileList[0].endswith(".yaml"):
            self.data_model.load_from_yaml(pyFileList[0])
        elif all( [f.endswith(".mzML") for f in pyFileList] ):
            self.data_model.loadFiles(pyFileList)
        else:
            # Check whether one of them is not mzML
            mzmls = [f for f in pyFileList if f.endswith("mzML")]
            others = [f for f in pyFileList if not f.endswith("mzML")]
            self.data_model.loadMixedFiles(mzmls, others)
        self._refresh_view(time=time.time()-start)

    def _refresh_view(self, time=0):
        """
        Refresh the whole application view (e.g. after loading new files)
        """

        # Get precursors from data model and then set the precursor tree structure
        precursor_model = self.application.get_precursor_model()
        precursor_model.set_precursor_tree_structure(self.data_model.get_precursor_tree())
        tmessage = ""
        if time > 0:
            tmessage = ". Loading took %0.4fs" % time
        self.statusBar().showMessage(self.data_model.getStatus() + tmessage)
        self.application.add_plots(self.data_model)

    def center(self):
        
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())
