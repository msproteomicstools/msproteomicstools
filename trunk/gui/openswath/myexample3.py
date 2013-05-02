
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
                QtGui.QColor( 255, 0, 255),
                QtGui.QColor( 255, 0, 255),
                QtGui.QColor( 255, 0, 255),
        ]
        self.curves = []

    def add_curves(self, nr):

        for i in range(nr):
            param = CurveParam()
            param.label = 'My curve'
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
        self.initialize()

    def initialize(self):
        self.colors = [
                QtGui.QColor( 255, 0, 0),
                QtGui.QColor( 50, 50, 50),
                QtGui.QColor( 255, 0, 255),
                QtGui.QColor( 255, 0, 255),
                QtGui.QColor( 255, 0, 255),
                QtGui.QColor( 255, 0, 255),
        ]
        self.curves = []
        self.ranges = []

    def add_curves(self, nr):

        plot = self.get_plot()
        for i in range(nr):
            param = CurveParam()
            param.label = 'My curve'
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

    def update_all_curves(self, data):
        for curve in self.curves:
            curve.set_data( range( 0, 20, 2), map( lambda _: random(), range( 0, 10 ) ) )
        # TODO only here we can get the range ???
        # print "range was ", self.myrange.get_range()
        r = int( random() * 15)
        self.myrange.set_range(r, r+5)

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
        
    def react_to_mouse(self):
        # print "react to mouse"
        pass

    def react_to_mouse_release(self):
        # print "react to release mouse"
        pass

    def initUI(self):
        self.layout = QtGui.QGridLayout(self)

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

class Communicate(QtCore.QObject):
    
    catch_mouse_press = QtCore.pyqtSignal() 
    catch_mouse_release = QtCore.pyqtSignal() 
    
class ApplicationView(QtGui.QWidget):
    
    def __init__(self):
        super(ApplicationView, self).__init__()
        
        self.initUI()
        
    def initUI(self):

        self.model = ExamplePeptides_()
        # self.model = ExamplePeptidesModel()
        self.model.setHorizontalHeaderLabels([self.tr("Peptides")])

        self.treeView = ExamplePeptidesTreeView()
        self.treeView.setModel(self.model)

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

        self.add_plots()

    def add_plots(self):
        
        self.plot = CurveDialogView(edit=False, toolbar=False )
        self.plot2 = CurvePlotView( self )

        self.plot.add_curves(3)
        self.plot2.add_curves(2)

        self.graph_layout.add_new(self.plot)
        self.graph_layout.add_new(self.plot2)

    def treeViewClicked(self, newvalue, oldvalue):
        print "got value", newvalue, "old", oldvalue
        self.plot.update_all_curves(None)
        self.plot2.update_all_curves(None)

    def widgetclicked(self, value):
        print "clicked iittt"


class MainWindow(QtGui.QMainWindow):
    
    def __init__(self):
        super(MainWindow, self).__init__()
        
        self.initUI()
        
    def initUI(self):               
        
        application = ApplicationView()
        self.setCentralWidget(application)

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

    def showDialog(self):

        fileList = QtGui.QFileDialog.getOpenFileNames(self, 'Open file')
        
        print "opened file ", fileList
        print "will open files", [str(f) for f in fileList]

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
