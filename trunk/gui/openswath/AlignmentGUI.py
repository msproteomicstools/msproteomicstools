#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
OpenSwath Viewer

Usage:
    middle-mouse click : auto-zoom
    middle-mouse move  : pan
    right-mouse move   : zoom

TODO: 
    use QDockWidget

INSTALL:
    Dependencies:
        - pymzML (from https://github.com/hroest/pymzML)
        - Python >= 2.5
        - numpy / scipy
        - PyQt4 >= 4.3
        - PyQwt >= 5.2
        - PIL (Python Imaging Library, maybe not necessary...)
        - guidata
        - guiqwt

    Install guiqwt and you should be fine, see https://code.google.com/p/guiqwt/

    - GNU/LINUX (Debian/Ubuntu/ArchLinux) should be packaged (python-guiqwt on ubuntu)
    - RedHat: PyQwt, PyQt, scipy, numpy should be the packages

    - Windows: install these dependencies:
                * http://pythonhosted.org/guiqwt/installation.html
                * maybe get started with http://pyqwt.sourceforge.net/download.html
"""

import sys,time

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, QModelIndex

from guiqwt.builder import make
from guiqwt.styles import CurveParam, COLORS
from guiqwt.transitional import QwtPlotItem

# global parameters
TITLE_FONT_SIZE = 10
AXIS_FONT_SIZE = 8
AUTOSCALE_Y_AXIS = True

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
from views.Plot import MultiLinePlot

# 
## The event handler if someone clicks on the graph 
#
class GraphEventHandler():
    """This object handles the events when the user zooms or pans.
    """
    def __init__(self, parent):
        self.plot = None
        self.parent = parent

    def panMouseMove(self, f, ev):
        pass

    def domove(self, f, ev):
        pass

    def panMousePress(self, f, ev):
        pass

    def panMouseRelease(self, f, ev):
        self.parent.reset_axis_all_plots(self.plot.get_plot().get_axis_limits("bottom"),
                                         self.plot.get_plot().get_axis_limits("left"), AUTOSCALE_Y_AXIS)

    def zoomMousePress(self, f, ev):
        pass
    def zoomMouseMove(self, f, ev):
        pass

    def zoomMouseRelease(self, f, ev):
        self.parent.reset_axis_all_plots(self.plot.get_plot().get_axis_limits("bottom"),
                                         self.plot.get_plot().get_axis_limits("left"), AUTOSCALE_Y_AXIS)

    def initialize(self, plot):
        self.plot = plot

        event_filter = self.plot.get_plot().filter
        try:
            # possible interactive choices: 
            # SelectTool is default, RectZoomTool needs toolbar, SignalStatsTool needs toolbar
            from guiqwt.tools import SelectTool
            [t for t in self.plot.tools
                if isinstance(t, SelectTool)]
            start_state = t[0].start_state.values()[0]
        except Exception:
            # fallback value (1 is the first state)
            start_state = 1
            # see the following call chain
            """
            CurveDialog.__init__
            -> CurveWidgetMixin.__init__
              -> register_tools
                -> register_all_curve_tools
                  -> register_standard_tools
                    -> add_tool
                      -> register_plot (e.g. SelectTool which is default)
                        -> setup_filter (e.g. SelectTool which is default)
                          -> new_state (StatefulEventFilter)
            """

        # For the Middle Button (pan)
        self.state0 = event_filter.add_event(start_state, event_filter.mouse_press(Qt.MidButton, Qt.NoModifier), self.panMousePress) 
        self.state1 = event_filter.add_event(self.state0, event_filter.mouse_move(Qt.MidButton, Qt.NoModifier), self.panMouseMove) 
        # event_filter.add_event(self.state1, event_filter.mouse_move(Qt.MidButton, Qt.NoModifier), self.domove) 
        event_filter.add_event(self.state1, event_filter.mouse_release(Qt.MidButton, Qt.NoModifier), self.panMouseRelease) 

        # For the Right Button (zoom)
        self.state10 = event_filter.add_event(start_state, event_filter.mouse_press(Qt.RightButton, Qt.NoModifier), self.zoomMousePress) 
        self.state11 = event_filter.add_event(self.state10, event_filter.mouse_move(Qt.RightButton, Qt.NoModifier), self.zoomMouseMove) 
        event_filter.add_event(self.state11, event_filter.mouse_release(Qt.RightButton, Qt.NoModifier), self.zoomMouseRelease) 

# 
## The widget for the Graphing area on the right
#
class GraphArea(QtGui.QWidget):

    def __init__(self):
        super(GraphArea, self).__init__()

        self.initUI()
        self._wcount = 1
        self.plots = []
        # self.c.catch_mouse_press.connect(self.react_to_mouse)
        # self.c.catch_mouse_release.connect(self.react_to_mouse_release)
         
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

        from guidata.qt.QtGui import QFont
        for i, run in enumerate(datamodel.get_runs()):

            self.plot = MultiLinePlot(edit=False, toolbar=False,
                                      options=dict(xlabel="Time (s)", ylabel="Intensity") )
            self.plot.setDataModel(run)

            # start event handler which will hook itself into the guiqwt event model
            e = GraphEventHandler(self)
            e.initialize(self.plot)

            # set font and title of plot
            self.plot.get_plot().font_title.setPointSize(TITLE_FONT_SIZE)
            self.plot.get_plot().set_title(run.get_id())
            ax_font = self.plot.get_plot().get_axis_font("left")
            ax_font.setPointSize(AXIS_FONT_SIZE)
            self.plot.get_plot().set_axis_font("left", ax_font)
            self.plot.get_plot().set_axis_font("bottom", ax_font)

            self.layout.addWidget(self.plot, i % 3, int(i/3) )
            self.plots.append(self.plot)

    def reset_axis_all_plots(self, x_range, y_range, autoscale_y=False):
        for i, pl in enumerate(self.plots):
            pl.set_x_limits(x_range[0], x_range[1])
            if autoscale_y:
                pl.set_y_limits_auto(x_range[0], x_range[1])
            else:
                pl.set_y_limits(y_range[0], y_range[1])
            pl.get_plot().replot()

    def update_all_plots(self, chr_transition):
        """
        We update the plots for all runs.
        """
        # Get all data, compute overall min/max
        xmins = []
        xmaxs = []
        pairs = []
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
            pl.update_all_curves(data, labels, ranges, mscore, intensity)
            pl.set_x_limits(min(xmins),max(xmaxs))
            pl.get_plot().replot()

#
## Main Widget
# 
class ApplicationView(QtGui.QWidget):
    
    def __init__(self, parent):
        super(ApplicationView, self).__init__()
        self.parent = parent
        self.initUI()
        
    def changeReturnPressedTest(self):
        print "return pressed with text", self.treeLineEdit.text()
        # TODO do something here!, e.g. cycle through the elements

    def changedTextTest(self, text):

        # TODO also allow the column to be set!
        # cmp source/VISUAL/SpectraViewWidget.C
        column = 2

        import re
        s = re.compile(str(text), re.IGNORECASE)
        m = self.treeView.model()
        for model_idx in self.treeView.iterTopLevelElements(column):
            display_data = m.data(model_idx, Qt.DisplayRole).toPyObject()
            if s.match(display_data):
                break
            
        if s.match(display_data):
            selectionModel = self.treeView.selectionModel()
            selectionModel.clearSelection()
            selectionModel.select(model_idx, QtGui.QItemSelectionModel.Select)
            self.treeView.setSelectionModel(selectionModel)
            # Now scroll to the item
            self.treeView.scrollTo(model_idx, QtGui.QAbstractItemView.PositionAtCenter)

    def initUI(self):

        self._precursor_model = PeptideTree([])
        self._precursor_model.setHorizontalHeaderLabels([self.tr("Peptides")])

        self.treeView = PeptidesTreeView()
        self.treeView.setModel(self._precursor_model)

        ## -> this would be how to sort but it seems to be buggy
        ## self.pProxyModel = QtGui.QSortFilterProxyModel()
        ## self.pProxyModel.setSourceModel(self._precursor_model)
        ## self.treeView.setModel(self.pProxyModel)
        ## self.treeView.setSortingEnabled(True)

        # Do the left side (hirarchical tree)
        self.leftside_layout = QtGui.QVBoxLayout()
        self.leftside_layout.addWidget(self.treeView)
        self.treeLineEdit = QtGui.QLineEdit()
        self.leftside_layout.addWidget(self.treeLineEdit)
        self.leftside = QtGui.QWidget()
        self.leftside.setLayout(self.leftside_layout)

        # TODO refactor to somewhere else!
        self.treeLineEdit.textChanged.connect(self.changedTextTest)
        self.treeLineEdit.returnPressed.connect(self.changeReturnPressedTest)

        # Do the main application (leftside/graphing area)
        self.graph_layout = GraphArea()
        horizontal_splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        horizontal_splitter.addWidget(self.leftside)
        horizontal_splitter.addWidget(self.graph_layout)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(horizontal_splitter)

        self.setLayout(hbox)

        # QItemSelectionModel -> connect the tree to here
        self.treeView.selectionModel().selectionChanged.connect(self.treeViewClicked) 

        # add dummy plots to the graph layout
        self.graph_layout.add_plots_dummy()

    def get_precursor_model(self):
        return self._precursor_model

    def set_communication(self, c):
        self.c = c

    def treeViewClicked(self, newvalue, oldvalue):
        if len(newvalue.indexes()) == 0 :
            return

        # assert that the the underlying selected element is always the same. 
        assert all(x.internalPointer() == newvalue.indexes()[0].internalPointer() for x in newvalue.indexes())

        # selected_precursor = newvalue.indexes()[0].internalPointer().ref.getName()
        import time
        s = time.time()
        self.graph_layout.update_all_plots(newvalue.indexes()[0].internalPointer().ref)
        self.parent.statusBar().showMessage(self.parent.data_model.getStatus() + ". Drawn plots in %0.4fs."  % (time.time() - s))


    def add_plots(self, datamodel):
        self.graph_layout.add_plots(datamodel)
        self.expandLevel("smart")

    def expandLevel(self, level):
        if level == "Peptides":
            self.treeView.expandToDepth(0)
        elif level == "Precursors":
            self.treeView.expandToDepth(1)
        elif level == "smart":
            self.treeView.expandMultiElementItems()
        else:
            self.treeView.collapseAll()

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

    def showDialog(self):

        fileList = QtGui.QFileDialog.getOpenFileNames(self, 'Open file')
        pyFileList = [str(f) for f in fileList]

        # Load the files
        start = time.time() 
        self.data_model.loadFiles(pyFileList)
        self._refresh_view(time=time.time()-start)

    def _refresh_view(self, time=0):

        # get precursors from data and set it 
        pr_list = self.data_model.get_precursor_list()
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
