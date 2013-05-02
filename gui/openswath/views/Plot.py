#!/usr/bin/python
# -*- coding: utf-8 -*-

from guiqwt.curve import CurveItem
from guiqwt.plot import CurvePlot, CurveDialog

from PyQt4 import QtGui 
from PyQt4.QtCore import Qt

from guiqwt.builder import make
from guiqwt.styles import CurveParam, COLORS
from guiqwt.transitional import QwtPlotItem

USE_ANTIALIASING = True

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

