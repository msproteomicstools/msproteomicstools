#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
OpenSwath Viewer
"""

import sys

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, QModelIndex

from guiqwt.plot import CurvePlot, CurveDialog
from guiqwt.curve import CurveItem
from guiqwt.builder import make
from guiqwt.styles import CurveParam, COLORS
from guiqwt.shapes import XRangeSelection

from guiqwt.transitional import QwtPlotItem

# global parameters
TITLE_FONT_SIZE = 10
AXIS_FONT_SIZE = 8
USE_ANTIALIASING = True
AUTOSCALE_Y_AXIS = True

class Communicate(QtCore.QObject):
    
    catch_mouse_press = QtCore.pyqtSignal() 
    catch_mouse_release = QtCore.pyqtSignal() 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class RunDataModel():

    def __init__(self, run, filename):
        import os
        self._run = run
        self._filename = filename
        self._basename = os.path.basename(filename)
        self._precursor_mapping = {}
        self._sequences_mapping = {}
        self._range_mapping = {}
        self._score_mapping = {}
        self._intensity_mapping = {}

        self._group_by_precursor()
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
                    # specific to pymzl, we need to get rid of those two entries
                    if key in ("indexList", "TIC"): continue

                    # The precursor identifier is the second (or third for
                    # decoys) element of the openswath format. Decoys should
                    # get a different identifier!
                    components = key.split("_")
                    trgr_nr = str(components[1])
                    if components[0].startswith("DECOY"):
                        trgr_nr = "DECOY_" + str(components[2])

                    if self._precursor_mapping.has_key(trgr_nr):
                        self._precursor_mapping[trgr_nr].append(key)
                    else:
                        self._precursor_mapping[trgr_nr] = [key]

        else:
            # TODO fallback option!!!
            pass
            raise Exception("Could not parse chromatogram ids ... ")

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

    def get_data_for_precursor(self, precursor):
        """Retrieve data for a specific precursor - data will be as list of
        pairs (timearray, intensityarray)"""

        if not self._precursor_mapping.has_key(str(precursor)):
            return [ [ [0], [0] ] ]

        transitions = []
        for chrom_id in self._precursor_mapping[str(precursor)]:
            c = self._run[str(chrom_id)] 
            transitions.append([c.time, c.i])

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
        return self._precursor_mapping.get(str(precursor), [])

    def get_precursors_for_sequence(self, sequence):
        return self._sequences_mapping.get(sequence, [])

    def get_all_precursor_ids(self):
        return self._precursor_mapping.keys()

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

class DataModel(object):

    def __init__(self):
        self.precursors = set([])
        self.runs = []

    def loadFiles(self, filenames):

        # load new files, clean up ...
        self.runs = []
        self.precursors = set([])

        for f in filenames:
            print "read file", f
            import pymzml
            run_ = pymzml.run.Reader(f, build_index_from_scratch=True)
            run = RunDataModel(run_, f)
            self.runs.append(run)
            self.precursors.update(run.get_all_precursor_ids())

    def getStatus(self):
        if len(self.runs) == 0:
            return "Ready"

        tr_cnt = 0
        for r in self.runs:
            tr_cnt += len(r._run.info['offsets']) -2

        return '%s Transitions, %s Peptides (total %s Transitions)' % ( 
            len(self.runs[0]._run.info['offsets']) -2, len(self.runs[0]._sequences_mapping), tr_cnt )

    def get_precursor_list(self):
        return self.precursors

    def get_precursor_tree(self):
        return self._build_tree()

    def _build_tree(self):

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

# 
## Tree Model - there are two options
#
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

    def getRange(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_range_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_range_data(prec[0]) 
        return [ 0,0]

    def getProbScore(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_score_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_score_data(prec[0]) 
        return 1.0

    def getIntensity(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_intensity_data(self.getName()) 
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_intensity_data(prec[0]) 
        return -1.0

    def getLabel(self, run):
        # return only the last element
        return [l.split("_")[-1] for l in self._getLabel(run)]

    def _getLabel(self, run):
        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_transitions_for_precursor(self.getName())
        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.name)
            if len(prec) == 1:
                return run.get_transitions_for_precursor(prec[0])
            else:
                pass
        elif CHROMTYPES[self.mytype] == "Transition" :
            return [self.getName()]
        return [ "" ]

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
        return 3

    def data(self, index, role):
        if not index.isValid():
            return None
        node = index.internalPointer()
        if role == Qt.DisplayRole and index.column() == 0:
            return QtCore.QVariant( node.ref.getPeptideSequence() )
        if role == Qt.DisplayRole and index.column() == 1:
            return QtCore.QVariant( node.ref.charge )
        if role == Qt.DisplayRole and index.column() == 2:
            return QtCore.QVariant( node.ref.name )
        return None

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 0:
            return 'Peptide Sequence'
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 1:
            return 'Charge'
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 2:
            return 'Full Name'
        return None

    def set_precursor_tree_structure(self, data):

        # first delete all rows
        parent = QModelIndex()
        self.beginRemoveRows(parent, 0, len(self.rootElements) )
        self.rootElements = []
        self.endRemoveRows()

        # now add the new rows
        parent = QModelIndex()
        self.beginInsertRows(parent, 0, len(data) )
        self.rootElements = data
        self.endInsertRows()

        # initialize super method again
        self.initialize()

    def setHorizontalHeaderLabels(self, parent):
        pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# 
## Tree View
#
class PeptidesTreeView( QtGui.QTreeView ):

    def __init__(self):
        super(PeptidesTreeView, self).__init__()

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        # self.customContextMenuRequested.connect(self.openMenu)
        
    def expandMultiElementItems(self):
        """
        Expand all top elements in the tree that have more than one child
        """
        for model_idx in self.iterTopLevelElements(0):
            if len(model_idx.internalPointer().subnodes) > 1:
                self.setExpanded(model_idx, True)

    def iterTopLevelElements(self, column):
        # find (fake) root and model
        root = self.rootIndex()
        m = self.model()
        for i in range(m.rowCount(root)):
            yield m.index(i,column,root)

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

    def create_curves(self, labels, this_range):

        self.curves = []
        plot = self.get_plot()
        plot.del_all_items(except_grid=False)
        for i,l in enumerate(labels):
            param = CurveParam()
            param.label = str(l)
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
            curve.setRenderHint(QwtPlotItem.RenderAntialiased, USE_ANTIALIASING)
            l = make.legend("TR")
            plot.add_item( l )

        self.myrange = make.range(this_range[0], this_range[1])
        self.myrange.itemChanged()
        self.myrange.set_movable(False)
        #
        if False:
            # make all colors black/gray
            # check /usr/lib/python2.7/dist-packages/guiqwt/config.py for parameters
            # 
            from guidata.dataset.dataitems import ColorItem
            self.myrange.shapeparam.fill = ColorItem("gray")
            self.myrange.shapeparam.line.color = "#000000"
            self.myrange.shapeparam.sel_line.color = "#000000"
            self.myrange.shapeparam.symbol.facecolor = "#000000"
            self.myrange.shapeparam.sel_symbol.facecolor = "#000000"
            #
            self.myrange.shapeparam.update_range(self.myrange) # creates all the above QObjects

        self.myrange.set_resizable(False)
        self.myrange.itemChanged()
        # print self.myrange._can_move, "move"
        # disp2 = make.computations(self.myrange, "TL",
        #                               [(curve, "min=%.5f", lambda x,y: y.min()),
        #                                (curve, "max=%.5f", lambda x,y: y.max()),
        #                                (curve, "avg=%.5f", lambda x,y: y.mean())])
        # plot.add_item( disp2 )
        plot.add_item( self.myrange )

    def rangeChanged(self, data):
        print "range changed"

    def set_x_limits(self, xmin, xmax):
        self.get_plot().set_axis_limits('bottom', xmin, xmax)

    def set_y_limits(self, ymin, ymax):
        self.get_plot().set_axis_limits('left', ymin, ymax)

    def set_y_limits_auto(self, xmin, xmax):
        """
        Automatically set y limits based on the data in the plot and the
        visible range xmin, xmax.
        """
        allmin = []
        allmax = []
        for curve in self.curves:
            data = curve.get_data()
            filtered_data = [y for x,y in zip(data[0], data[1]) if x>xmin and x<xmax]
            if len(filtered_data) == 0:
                continue
            allmin.append( min(filtered_data) )
            allmax.append( max(filtered_data) )
        if len(allmin) == 0 or len(allmax) == 0 : return
        self.set_y_limits(min(allmin), max(allmax) )

    def update_all_curves(self, data, labels, ranges, mscore, intensity):

        assert len(data) == len(labels)
        self.create_curves(labels, ranges)
        if mscore is not None: 
            self.mscore_label = make.label("m_score=%0.4g" % mscore, "TL", (0,0), "TL")
            self.get_plot().add_item( self.mscore_label )
            self.l2 = make.label("Int=%0.4g" % intensity, "TL", (0,25), "TL")
            self.get_plot().add_item( self.l2 )
            self.width = make.label("PeakWidth=%0.4g" % (ranges[1]-ranges[0]), "TL", (0,50), "TL")
            self.get_plot().add_item( self.width )

        for d, curve in zip(data, self.curves):
            curve.set_data( d[0], d[1] )

        self.get_plot().do_autoscale()

    def mouseReleaseEvent(self, event):
        pass
        # TODO capture the mouse release 
        # print "mouse was released, range was ", self.myrange.get_range()

class GraphEventHandler():
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
        self.data_model.loadFiles(pyFileList)
        self._refresh_view()

    def _refresh_view(self):

        # get precursors from data and set it 
        pr_list = self.data_model.get_precursor_list()
        precursor_model = self.application.get_precursor_model()
        precursor_model.set_precursor_tree_structure(self.data_model.get_precursor_tree())
        self.statusBar().showMessage(self.data_model.getStatus())
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
