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

import sys,time, re
import argparse

## We need to import guiqt first
## https://github.com/pierreraybaut/guidata/issues/35#issuecomment-171573060
have_guiqwt = True
try:
    from guiqwt.curve import CurveItem
except ImportError:
    print "Could not import guiqwt, will try to use Qwt only."
    have_guiqwt = False

# from guidata import qt
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, QModelIndex

# global parameters
TITLE_FONT_SIZE = 10
AXIS_FONT_SIZE = 8
# Turning off guiqwt (USE_GUIQWT=False) should be safer, it then uses the
# fallback to plain Qwt.
USE_GUIQWT = False
USE_GUIQWT = True

class Communicate(QtCore.QObject):
    
    catch_mouse_press = QtCore.pyqtSignal() 
    catch_mouse_release = QtCore.pyqtSignal() 

# 
## Data Model and Tree Model
#
from openswathgui import *
from openswathgui.models.MSData import DataModel
from openswathgui.models.PeptideTree import PeptideTree

# 
## Views for the plots and the peptide tree (on the right)
#
from openswathgui.views.PeptideTree import PeptidesTreeView
from openswathgui.views.Plot import have_guiqwt 

# 
## The widget for the graphing area on the right
#
class GraphArea(QtGui.QWidget):
    """
    The Graph Area is displayed on the right side of the main area (see :class:`.ApplicationView`)

    The actual implementation of the plotting is in the Plot model, see :mod:`.openswathgui.views.Plot`

    Attributes:
        self.plots: The underlying list of plots (either of type :class:`.GuiQwtMultiLinePlot` or :class:`.QwtMultiLinePlot`)
    """

    def __init__(self, use_guiqwt = USE_GUIQWT):
        super(GraphArea, self).__init__()

        self._initUI()
        self._wcount = 1
        self.nr_rows = 3
        self.autoscale_y_axis = True
        self.plots = []

        self.changePlotEngine(use_guiqwt)

        # self.c.catch_mouse_press.connect(self.react_to_mouse)
        # self.c.catch_mouse_release.connect(self.react_to_mouse_release)

    def changePlotEngine(self, use_guiqwt):
        if use_guiqwt and have_guiqwt:
            from openswathgui.views.Plot import GuiQwtMultiLinePlot as MultiLinePlot
            self.MultiLinePlot = MultiLinePlot
        else:
            from openswathgui.views.Plot import QwtMultiLinePlot as MultiLinePlot
            self.MultiLinePlot = MultiLinePlot
         
    @QtCore.pyqtSlot(float, float, float, float)
    def plotZoomChanged(self, xmin, xmax, ymin, ymax):
        """
        Slot to deal with the underlying plot emitting a zoomChanged signal
        """
        # after (moving the image), adjust and replot _all_ plots
        self._reset_axis_all_plots([xmin, xmax], [ymin, ymax], self.autoscale_y_axis)

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

    def _initUI(self):
        self.layout = QtGui.QGridLayout(self)

    def _delete_all(self):
        for i in range(self.layout.count()):
            self.layout.itemAt(i).widget().close()
        self._wcount = 1

    def _add_new(self, l):
        self.layout.addWidget(l, self._wcount, 0)
        self._wcount += 1

    def add_plots_dummy(self):
        """
        Add dummy plots for testing
        """
        
        self.plots = []

        self.plot = self.MultiLinePlot(edit=False, toolbar=False )
        self.plot.create_curves([1,2,3], [ [0,0] ] )
        self._add_new(self.plot)
        self.plots.append(self.plot)

        #self.plot2 = CurvePlotView( self )
        self.plot2 = self.MultiLinePlot(edit=False, toolbar=False )
        self.plot2.create_curves([1,2], [ [0,0] ])
        self._add_new(self.plot2)
        self.plots.append(self.plot2)

    def add_plots(self, datamodel):
        """
        Add a plot for each run that needs to be displayed

        Args:
            datamodel(:class:`.DataModel`): The data model containing data to plot
        """
        
        self.plots = []
        self._delete_all()

        for i, run in enumerate(datamodel.get_runs()):

            self.plot = self.MultiLinePlot(edit=False, toolbar=False,
                                      options=dict(xlabel="Time (s)", ylabel="Intensity") )
            self.plot.setDataModel(run)

            # set font and title of plot
            self.plot.setTitleFontSize(TITLE_FONT_SIZE)
            self.plot.setAxisFontSize(AXIS_FONT_SIZE)

            self.plot.zoomChanged.connect(self.plotZoomChanged) 

            self.layout.addWidget(self.plot, i % self.nr_rows, int(i/self.nr_rows) )
            self.plots.append(self.plot)

    def _reset_axis_all_plots(self, x_range, y_range, autoscale_y_axis=False):
        for i, pl in enumerate(self.plots):
            pl.set_x_limits(x_range[0], x_range[1])
            if autoscale_y_axis:
                pl.set_y_limits_auto(x_range[0], x_range[1])
            else:
                pl.set_y_limits(y_range[0], y_range[1])
            pl.replot()

    def update_all_plots(self, chr_transition, show_legend):
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
            assay_rt = chr_transition.getAssayRT(pl.run) 

            xminsr = [r[0] for r in ranges]
            xmaxsr = [r[1] for r in ranges]
            xwidth = [r[1] - r[0] for r in ranges]

            # print "Got data for", i, "plot labels", labels, "from", chr_transition, "and nr data ", len(data)
            # print "Got mscore and int", mscore, intensity
            # print "Got ranges ", ranges, labels, mscore, intensity

            # this next command takes about 10 ms per plot with Qwt, ca 30-40 ms with GuiQwt
            pl.update_all_curves(data, labels, ranges, mscore, intensity, assay_rt, show_legend)

            if mscore is None:
                pl.set_x_limits(min(xmins),max(xmaxs))
            else:
                pl.set_x_limits(min(xminsr) - max(xwidth),max(xmaxsr) + max(xwidth))
                pl.set_y_limits_auto(min(xminsr) - max(xwidth),max(xmaxsr) + max(xwidth))

            pl.replot()

#
## Peptide Tree Widget (left side)
# 
class PeptideTreeWidget(QtGui.QWidget):
    """
    The Peptide Tree Widget is displayed on the left side of the main area (see
    :class:`.ApplicationView`). It consists of a :class:`.PeptidesTreeView`
    widget and a search box below.

    Attributes:
        - self._precursor_model: The underlying peptide tree model (of type :class:`.PeptideTree`)
        - self.treeView: The underlying peptide tree view widget (of type :class:`.PeptidesTreeView`)
        - self.treeLineEdit: A place for input below the tree
        - self.treeComboBox: A combo box relating to the column for search

    Emits the following signals:
        - selectionChanged : when the peptide selection is changed
    """

    # Signals
    selectionChanged = QtCore.pyqtSignal(QModelIndex)
    """
    Qt signal emitted when the peptide selection changes
    """

    def __init__(self, firstColumnName):
        super(PeptideTreeWidget, self).__init__()

        self._precursor_model = None
        self.treeView = None

        self.first_column_name_ = firstColumnName
        self.initUI()

    def initUI(self):
        """
        Set up the model and the view

        This sets up the layout as follows:
            1. self.treeView widget
            2. searchbox_layout layout
                2.1 self.treeLineEdit (QtGui.QLineEdit)
                2.2 self.treeComboBox (QtGui.QComboBox)
        """

        self._precursor_model = PeptideTree([], firstColumnName=self.first_column_name_)
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
        final_layout = QtGui.QVBoxLayout()
        final_layout.addWidget(self.treeView)
        final_layout.addLayout(searchbox_layout)
        self.setLayout(final_layout)

        #
        ## Connect the Signals
        #
        self.treeLineEdit.textChanged.connect(self.changedText)
        self.treeLineEdit.returnPressed.connect(self.changeReturnPressedSlot)
        # QItemSelectionModel -> connect the tree to here
        self.treeView.selectionModel().selectionChanged.connect(self.treeViewSelectionChangedSlot) 

    def changeReturnPressedSlot(self):
        """
        Slot connected to the signal generated by pressing return in the text field
        """

        if self.treeiter is None:
            return

        column =  self.treeComboBox.currentIndex()
        try:
            model_idx = self.treeiter.next()
        except StopIteration:
            # If we have reached the bottom, wrap around
            try:
                self.treeiter = self._iterIndices(self.treeView, column, self.treeLineEdit.text())
                model_idx = self.treeiter.next()
            except StopIteration:
                # No match found
                return

        self.treeView.selectAndScrollTo(model_idx)

    def _iterIndices(self, treeView, column_, text_):
        """
        Iterator generator function to go through all indexes
        """

        s = re.compile(str(text_), re.IGNORECASE)
        m = treeView.model()

        for model_idx in treeView.iterAllLevelElements(column_):
            display_data = m.data(model_idx, Qt.DisplayRole)
            if s.search(display_data):
                yield model_idx

    def changedText(self, text):
        """
        Slot connected to the signal generated by changing the text in the text field
        """

        column =  self.treeComboBox.currentIndex()
        try:
            self.treeiter = self._iterIndices(self.treeView, column, text)
            model_idx = self.treeiter.next()
        except StopIteration:
            # No match found
            return

        self.treeView.selectAndScrollTo(model_idx)

    @QtCore.pyqtSlot(QtGui.QItemSelectionModel, QtGui.QItemSelectionModel)
    def treeViewSelectionChangedSlot(self, newvalue, oldvalue):
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
        """
        Access to the underlying precursor model

        Returns
        -------
        precursor_model : :class:`.PeptideTree`
            The underlying precursor model 
        """
        return self._precursor_model


#
## Main Widget
# 
class ApplicationView(QtGui.QWidget):

    """
    The main/central widget for the application which is directly called from the MainWindow

    Attributes:
        self.leftside: Reference to the left side widget (of type :class:`.PeptideTreeWidget`)
        self.graph_layout: Reference to the right side widget (of type :class:`.GraphArea`)

    Emits the following signals:
        - plotsUpdate : when the plots need to be updated
    """
    
    # Signals
    plotsUpdated = QtCore.pyqtSignal(float)
    """
    Qt signal emitted when plots need to be updated due to selecting a different peptide
    """

    def __init__(self, parent, settings):
        super(ApplicationView, self).__init__()
        self.parent = None
        self.treeiter = None
        self.settings = settings

        self.leftside = None # leftside widget
        self.graph_layout = None # rightside widget (graph area)

        self._initUI()
        
    @QtCore.pyqtSlot(QModelIndex)
    def treeSelectionChanged(self, idx):
        """
        Grab the selectionChanged signal from the :class:`.PeptideTreeWidget` and
        accordingly update the graphing area on the right.
        """
        s = time.time()
        self.graph_layout.update_all_plots(idx.internalPointer().ref, self.settings.show_legend)
        self.plotsUpdated.emit(time.time()-s)

    def _initUI(self):

        # Do the peptide tree on the left side and connect its signals
        self.leftside = PeptideTreeWidget(self.settings.first_column_name_)
        self.leftside.selectionChanged.connect(self.treeSelectionChanged)

        # Do the main application (graphing area on the right side)
        self.graph_layout = GraphArea(settings.use_guiqwt)
        horizontal_splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        horizontal_splitter.addWidget(self.leftside)
        horizontal_splitter.addWidget(self.graph_layout)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(horizontal_splitter)

        self.setLayout(hbox)

        # add dummy plots to the graph layout
        self.graph_layout.add_plots_dummy()

    def get_precursor_model(self):
        """
        Access to the underlying precursor model

        Returns
        -------
        precursor_model : The underlying precursor model of class :class:`.PeptideTree`
        """
        return self.leftside.get_precursor_model()

    def set_communication(self, c):
        self.c = c

    def add_plots(self, datamodel):
        """
        Add a plot for each run that needs to be displayed (calls the underlying :class:`.GraphArea`)

        Args:
            datamodel(:class:`.DataModel`): The data model containing data to plot
        """
        self.graph_layout.add_plots(datamodel)
        self.leftside.expandLevel("smart")

    def widgetclicked(self, value):
        print "clicked iittt"


#
## Settings object
# 
class Settings(object):

    def __init__(self, runMode, use_guiqwt = False):
        self.show_legend = True
        self.draw_transitions = False
        self.autoscale_y_axis = True
        self.nr_rows = 3
        self.window_title = 'TAPIR'
        self.use_guiqwt = use_guiqwt

        if runMode == "proteomics":
            self.first_column_name_ = 'Identifier'
        elif runMode == "metabolomics":
            self.first_column_name_ = 'Metabolite'
        else:
            raise Exception("Unknown run mode %s" % runMode)

#
## Configuration Dialog
# 
class ConfigDialog(QtGui.QDialog):
    """
    Configuration Dialog 
    """

    def __init__(self, parent, settings):
        QtGui.QDialog.__init__(self, parent)
        self.settings = settings
        self.parent = parent
        self._initUI()

    def closeAndSave(self):
        """
        Responds to the close and save action
        """

        if not have_guiqwt and self.use_guiqwt.isChecked():
            box = QtGui.QMessageBox.about(self, "Error", "QtGui is not available, please remove it.")
            box.setIcont(QtGui.QMessageBox.Critical)
            return

        self.settings.window_title = str(self.window_title.text())
        self.settings.nr_rows = int(self.nr_rows.text())
        self.settings.show_legend = self.show_legend.isChecked()
        self.settings.draw_transitions = self.draw_transitions.isChecked()
        self.settings.autoscale_y_axis = self.autoscale_y_axis.isChecked()
        self.settings.use_guiqwt = self.use_guiqwt.isChecked()
        self.parent.updateSettings(self.settings)
        self.close()

    def _initUI(self):

        # Right side layout
        # contentsWidget = QtGui.QListWidget()

        # Close button
        self.closeButton = QtGui.QPushButton("Close");
        self.closeButton.clicked.connect(self.closeAndSave)
        # connect(closeButton, SIGNAL(clicked()), this, SLOT(close()));

        # Left side layout
        updateGroup = QtGui.QGroupBox("TAPIR settings");
        self.show_legend = QtGui.QCheckBox("Show legend");
        self.draw_transitions = QtGui.QCheckBox("Draw individual transitions");
        self.autoscale_y_axis = QtGui.QCheckBox("Autoscale y axis");
        self.use_guiqwt = QtGui.QCheckBox("Use guiqwt advanced graphics (slower)");
        label_rows = QtGui.QLabel("Number of window rows");
        self.nr_rows = QtGui.QLineEdit();
        label_rows = QtGui.QLabel("Window Title");
        self.window_title = QtGui.QLineEdit();

        self.show_legend.setChecked( self.settings.show_legend )
        self.draw_transitions.setChecked( self.settings.draw_transitions )
        self.autoscale_y_axis.setChecked( self.settings.autoscale_y_axis )
        self.use_guiqwt.setChecked( self.settings.use_guiqwt )
        self.nr_rows.setText( str(self.settings.nr_rows) )
        self.window_title.setText( str(self.settings.window_title) )

        updateLayout = QtGui.QVBoxLayout()
        updateLayout.addWidget(self.show_legend);
        updateLayout.addWidget(self.draw_transitions);
        updateLayout.addWidget(self.autoscale_y_axis);
        updateLayout.addWidget(self.use_guiqwt);
        updateLayout.addWidget(label_rows);
        updateLayout.addWidget(self.nr_rows);
        updateLayout.addWidget(self.window_title);
        updateGroup.setLayout(updateLayout);

        ###################################
        # Composition of elements
        ###################################
        self.horizontalLayout = QtGui.QHBoxLayout()
        # self.horizontalLayout.addWidget(contentsWidget);
        self.horizontalLayout.addWidget(updateGroup);

        self.buttonsLayout = QtGui.QHBoxLayout()
        self.buttonsLayout.addStretch(1);
        self.buttonsLayout.addWidget(self.closeButton);

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addLayout(self.horizontalLayout);
        mainLayout.addStretch(1);
        mainLayout.addSpacing(12);
        mainLayout.addLayout(self.buttonsLayout);
        self.setLayout(mainLayout);

        self.setWindowTitle("Config Dialog");

#
## Main Window
# 
class MainWindow(QtGui.QMainWindow):
    """
    The main window running the application.

    - It contains a reference to the actual MS data model (self.data_model [models/MSData.py])
    - It contains a reference to the main widget (self.application_view [ApplicationView])

    - It loads files through self.loadFiles which delegates the call to the data model

    Attributes:
        self.data_model: Reference to the underlying data model (of type :class:`.DataModel`)
        self.application: Reference to the actual main widget (of type :class:`.ApplicationView`)
    """
    
    def __init__(self, settings):
        super(MainWindow, self).__init__()
        
        self.settings = settings
        self.c = Communicate()
        self.data_model = DataModel() # see models/MSData.py
        self.application = None # main application (ApplicationView), see initUI

        self._initUI()
        
    def _initUI(self):               

        self.data_model.setDrawTransitions( self.settings.draw_transitions )
        
        self.application = ApplicationView(self, self.settings)
        self.application.set_communication(self.c)
        self.application.plotsUpdated.connect(self.plotsUpdated)
        self.setCentralWidget(self.application)

        ###################################
        # Actions
        ###################################
        openFileIcon = QtGui.QIcon()
        openFileIcon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_DirOpenIcon),
                QtGui.QIcon.Normal, QtGui.QIcon.On)
        openFile = QtGui.QAction(openFileIcon, 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.showFileLoadDialog)

        openSettingsIcon = QtGui.QIcon("")
        openSettingsIcon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_FileIcon),
                QtGui.QIcon.Normal, QtGui.QIcon.On)
        openSettings = QtGui.QAction(openSettingsIcon, 'Open Settings', self)
        openSettings.setStatusTip('Open settings dialog')
        openSettings.triggered.connect(self._showSettings)

        exitIcon = QtGui.QIcon("")
        exitAction = QtGui.QAction(exitIcon, 'Exit', self)
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
        fileMenu.addAction(openSettings)       
        fileMenu.addAction(exitAction)

        ###################################
        # Toolbar
        ###################################
        toolbar = self.addToolBar('Exit')
        toolbar.addAction(openFile)
        toolbar.addAction(openSettings)
        # toolbar.addAction(exitAction)
        
        # self.setGeometry(300, 300, 250, 150)
        self.resize(850, 550)
        self._center()
        self.setWindowTitle(self.settings.window_title)
        self.show()
        self.statusBar().showMessage('Ready')

    @QtCore.pyqtSlot(float)
    def plotsUpdated(self, time_taken):
        """
        Qt slot: updates the status bar when :class:`.ApplicationView` emits :meth:`.ApplicationView.plotsUpdated`
        """
        self.statusBar().showMessage(self.data_model.getStatus() + ". Drawn plots in %0.4fs."  % (time_taken))

    def showFileLoadDialog(self):
        """
        Show the open file dialog to load a new dataset.
        """

        fileList = QtGui.QFileDialog.getOpenFileNames(self, 'Open dataset')
        self.loadFiles( [str(f) for f in fileList] )

    def loadFiles(self, pyFileList, fileType=None):
        """ Load a set of files to display

        1. Try to load single yaml file
        2. Try to load a list of only mzML files
        3. Try to load a mixed list of mzML and other files (.tsv)

        For the third option, 

        pyFileList : list of str
            List of paths to files
        fileType : str
            Description of the type of file the metadata file (valid: simple, yaml, traml, openswath)
        """

        start = time.time() 

        if fileType == "yaml" and len(pyFileList) == 1:
            Exception("When providing a yaml file, please do not provide any other files as input.")

        if len(pyFileList) == 1 and (pyFileList[0].endswith(".yaml") or fileType == "yaml"):
            self.data_model.load_from_yaml(pyFileList[0])
        elif all( [f.lower().endswith("mzml") for f in pyFileList] ):
            self.data_model.loadFiles(pyFileList)
        elif all( [f.lower().endswith("sqmass") for f in pyFileList] ):
            self.data_model.loadSqMassFiles(pyFileList)
        else:
            if any( [f.lower().endswith("sqmass") for f in pyFileList] ):

                # Separate the mzML and other files
                sqmass = [f for f in pyFileList if f.lower().endswith("sqmass")]
                others = [f for f in pyFileList if not f.lower().endswith("sqmass")]

                fileType = "sqmass"
                self.data_model.loadMixedFiles(sqmass, others, fileType)
            else:

                # Separate the mzML and other files
                mzmls = [f for f in pyFileList if f.lower().endswith("mzml")]
                others = [f for f in pyFileList if not f.lower().endswith("mzml")]

                if fileType == None and len(others) == 1 and others[0].lower().endswith("traml"):
                    fileType = "traml"

                self.data_model.loadMixedFiles(mzmls, others, fileType)

        # After loading, refresh view and print load time
        self._refresh_view(time=time.time()-start)

    def updateSettings(self, settings):
        """
        Update global settings (after closing the settings dialog)
        """

        self.application.graph_layout.nr_rows = settings.nr_rows
        self.application.graph_layout.autoscale_y_axis = settings.autoscale_y_axis
        self.application.graph_layout.changePlotEngine(settings.use_guiqwt)
        self.data_model.setDrawTransitions( settings.draw_transitions )
        self.setWindowTitle(self.settings.window_title)
        self._refresh_view()

    def _showSettings(self):
        settings = ConfigDialog(self, self.settings)
        settings.show()

    def _refresh_view(self, time=0):
        """
        Refresh the whole application view (e.g. after loading new files or changing the view)
        """

        # Get precursors from data model and then set the precursor tree structure
        precursor_model = self.application.get_precursor_model()
        precursor_model.set_precursor_tree_structure(self.data_model.get_precursor_tree())
        tmessage = ""
        if time > 0:
            tmessage = ". Loading took %0.4fs" % time

        self.statusBar().showMessage(self.data_model.getStatus() + tmessage)
        self.application.add_plots(self.data_model)

    def _center(self):
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

def handle_args():
    usage = ""
    usage += "\nThis program will display chromatograms and associated peakgroups."

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", required=False, nargs = '+', 
                        help = 'A list of input files (.chrom.mzML and feature_alignment output files).')
    parser.add_argument('--fileType', dest="filetype", required=False, 
                        help = 'Type of files describing the relations (simple, openswath, yaml)')
    parser.add_argument('--runMode', dest="run_mode", required=False, default="proteomics",
                        help = 'Mode to run in (proteomics, metabolomics)')
    parser.add_argument('--use_guiqwt', dest="use_guiqwt", required=False, default="False",
                        help = 'Whether to use guiQwt (False,True)')
    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':

    # Handle command line options
    options = handle_args()
    settings = Settings(options.run_mode, options.use_guiqwt == "True")

    # Set up Qt application
    app = QtGui.QApplication(sys.argv)
    ex = MainWindow(settings)

    # Check whether any options were given on the commandline
    if options.infiles is not None:
        ex.loadFiles(options.infiles, options.filetype)

    sys.exit(app.exec_())

