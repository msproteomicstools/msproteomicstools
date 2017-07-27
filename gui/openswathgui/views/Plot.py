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

from PyQt4 import QtGui, Qt, QtCore
from PyQt4.Qwt5 import QwtPlotItem
import PyQt4.Qwt5 as Qwt

# There are two implementations of the plotting view, one uses guiqwt (may not
# be present on all systems) and the other one uses plain Qwt (should be safer)
have_guiqwt = True
try:
    from guiqwt.curve import CurveItem
    from guiqwt.plot import CurvePlot, CurveDialog
    from guiqwt.builder import make
    from guiqwt.styles import CurveParam, COLORS
    from guiqwt.transitional import QwtPlotItem
except ImportError:
    print "Could not import guiqwt, will try to use Qwt only."
    have_guiqwt = False
    class CurveItem(object): pass
    class CurveDialog(object): 
        def __init__(self, *args, **kwargs):
            pass

try:
    QtGui.QStaticText("Test")
except AttributeError:
    print "Could not find QtGui.QStaticText - maybe your version of Qt and PyQt is too old? Will not be able to display text on top of the graphs."

USE_ANTIALIASING = True

class CurveItemModel(CurveItem):
    """
    A single curve
    """

    def __init__(self, *args, **kwargs):
        super(CurveItemModel, self).__init__(*args, **kwargs)

# 
## The widget for a single plot on the right
#
class GuiQwtMultiLinePlot(CurveDialog):
    """
    Widget for an individual plot of transitions

    For the Curve window we could use a CurveDialog or a CurvePlot. 

    CurveDialog has more features and seems more advanced.
    """

    # Signals
    zoomChanged = QtCore.pyqtSignal(float, float, float, float)

    def __init__(self, *args, **kwargs):
        super(GuiQwtMultiLinePlot, self).__init__(*args, **kwargs)
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

    def setTitleFontSize(self, fontsize):
        self.get_plot().font_title.setPointSize(fontsize)
        self.get_plot().set_title(self.run.get_id())

    def setAxisFontSize(self, fontsize):
        ax_font = self.get_plot().get_axis_font("left")
        ax_font.setPointSize(fontsize)
        self.get_plot().set_axis_font("left", ax_font)
        self.get_plot().set_axis_font("bottom", ax_font)

    def create_curves(self, labels, this_range, show_legend=True):

        self.curves = []
        plot = self.get_plot()
        plot.del_all_items(except_grid=False)
        for i,l in enumerate(labels):
            param = CurveParam()
            param.label = str(l)
            color = COLORS.get(self.colors[i % len(self.colors)],  self.colors[i % len(self.colors)] )
            param.line.color = color

            # create a new curve
            curve = CurveItemModel(param)
            self.curves.append(curve)
            plot.add_item( curve )
            curve.setRenderHint(QwtPlotItem.RenderAntialiased, USE_ANTIALIASING)
            l = make.legend("TR")
            if show_legend:
                plot.add_item( l )

        self.myranges = []
        for r in this_range:
            self.add_range_to_plot(plot, r)

    def add_range_to_plot(self, plot, this_range):

        self.myrange = make.range(this_range[0], this_range[1])
        self.myrange.itemChanged()
        self.myrange.set_movable(False)

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
        self.myranges.append( self.myrange )

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

    def update_all_curves(self, data, labels, ranges, mscore, intensity, assay_rt, show_legend):

        assert len(data) == len(labels), "Length of labels %s does not correspond to length of data %s" % (len(labels), len(data))
        self.create_curves(labels, ranges, show_legend)

        if mscore is not None and len(ranges) < 2: 
            self.mscore_label = make.label("m_score=%0.4g" % mscore, "TL", (0,0), "TL")
            self.get_plot().add_item( self.mscore_label )
            self.l2 = make.label("Int=%0.4g" % intensity, "TL", (0,25), "TL")
            self.get_plot().add_item( self.l2 )

            self.width = make.label("PeakWidth=%0.3fs" % (ranges[0][1]-ranges[0][0]), "TL", (0,50), "TL")
            self.get_plot().add_item( self.width )
            self.assay_rt = make.label("Assay RT=%0.3fs" % (assay_rt), "TL", (0,75), "TL")
            self.get_plot().add_item( self.assay_rt )

        for d, curve in zip(data, self.curves):
            curve.set_data( d[0], d[1] )

        self.get_plot().do_autoscale()

    def mouseReleaseEvent(self, event):
        pass

    def replot(self):
        self.get_plot().replot()

    def mouseReleaseEvent(self, event):
        xaxis_limits = self.get_plot().get_axis_limits("bottom")
        yaxis_limits = self.get_plot().get_axis_limits("left")
        self.zoomChanged.emit(xaxis_limits[0], xaxis_limits[1], yaxis_limits[0], yaxis_limits[1])

# 
## The widget for a single plot on the right using Qwt only
#
class QwtMultiLinePlot(Qwt.QwtPlot):
    """Widget for an individual plot of transitions (using Qwt Plot)

    Implements the same interface as the GuiQwtMultiLinePlot, but performs faster.
    """

    # Signals
    zoomChanged = QtCore.pyqtSignal(float, float, float, float)

    def __init__(self, *args, **kwargs):
        super(QwtMultiLinePlot, self).__init__(*args)
        self.myrange = None
        self.run = None
        self.initialize()
        self.has_mscore = False
        self.labels = []

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

        self.setCanvasBackground(Qt.Qt.white)

        picker_on = Qwt.QwtPicker.AlwaysOn
        picker_on = Qwt.QwtPicker.AlwaysOff
        self.zoomer = Qwt.QwtPlotZoomer(Qwt.QwtPlot.xBottom,
                                   Qwt.QwtPlot.yLeft,
                                   Qwt.QwtPicker.DragSelection,
                                   picker_on,
                                   self.canvas())
        self.zoomer.setRubberBandPen(QtGui.QPen(Qt.Qt.black))
        self.zoomer.setTrackerPen(QtGui.QPen(Qt.Qt.black))

        self.panner = Qwt.QwtPlotPanner( self.canvas() )
        self.panner.setMouseButton(Qt.Qt.MidButton)

    def create_curves(self, labels, this_range, show_legend=True):

        # delete / detach all old curves
        for c in self.curves:
            c.detach()
            
        # create new curves
        self.curves = []
        for i,curve_label in enumerate(labels):

            # Create legend item and set font
            curve_label_t = Qwt.QwtText(str(curve_label))
            f = curve_label_t.font()
            f.setPointSize(8)
            curve_label_t.setFont(f)

            curve = Qwt.QwtPlotCurve(curve_label_t)
            curve.attach(self)
            curve.setPen(Qt.QPen(self.colors[i % len(self.colors)]))
            curve.setRenderHint(QwtPlotItem.RenderAntialiased, USE_ANTIALIASING)
            self.curves.append(curve)

        if show_legend:
            legend = Qwt.QwtLegend()
            self.insertLegend(legend, Qwt.QwtPlot.BottomLegend);
            # self.insertLegend(legend, Qwt.QwtPlot.ExternalLegend);
        xaxis_title = Qwt.QwtText("Time (seconds)")
        yaxis_title = Qwt.QwtText("Intensity")
        self.setAxisTitle(Qwt.QwtPlot.xBottom, xaxis_title)
        self.setAxisTitle(Qwt.QwtPlot.yLeft, yaxis_title)

    def clearZoomStack(self):
        """Auto scale and clear the zoom stack
        """
        self.setAxisAutoScale(Qwt.QwtPlot.xBottom)
        self.setAxisAutoScale(Qwt.QwtPlot.yLeft)
        self.replot()
        self.zoomer.setZoomBase()

    def setDataModel(self, run):
        self.run = run

    def setTitleFontSize(self, fontsize):
        title = Qwt.QwtText(self.run.get_id() )
        titlefont = title.font()
        titlefont.setPointSize(fontsize)
        title.setFont(titlefont)
        self.setTitle(title)

    def setAxisFontSize(self, fontsize):
        ax_font = self.axisFont(Qwt.QwtPlot.xBottom)
        ax_font.setPointSize(fontsize)
        self.setAxisFont(Qwt.QwtPlot.xBottom, ax_font)
        self.setAxisFont(Qwt.QwtPlot.yLeft, ax_font)
        self.replot()

    def drawCanvas(self, painter):
        super(QwtMultiLinePlot, self).drawCanvas(painter)
        if not self.has_mscore: 
            return

        # add the labels
        current_height = 0
        for label_txt in self.labels:
            painter.drawStaticText ( 10, 10 + current_height, label_txt);
            current_height +=  label_txt.size().height()

        # add the range
        self.draw_range(painter, self.l_width, self.r_width)

    def draw_range(self, painter, l_width, r_width): 
        #, xMap, yMap, canvasRect):
        xMap = self.canvasMap(Qwt.QwtPlot.xBottom)
        yMap = self.canvasMap(Qwt.QwtPlot.yLeft)
        rct = self.canvas().contentsRect()

        rct2 = QtCore.QRectF(rct)
        rct2.setLeft(xMap.transform(l_width))
        rct2.setRight(xMap.transform(r_width))

        # TODO abstract this
        transparency = 40 # alpha channel  / transparency
        col =  QtGui.QColor( 255, 0, 0, transparency)   
        self.brush = QtGui.QBrush(col)
        pen = QtGui.QPen(col)

        # paint the filling rectable and the two lines left/right
        painter.fillRect(rct2, self.brush)
        painter.setPen(pen)
        painter.drawLine(rct2.topRight(), rct2.bottomRight())
        painter.drawLine(rct2.topLeft(), rct2.bottomLeft())

        # draw dashed line in the center
        dash = QtGui.QPen(pen)
        dash.setStyle(QtCore.Qt.DashLine)
        dash.setWidth(1)
        painter.setPen(dash)
        painter.drawLine(rct2.center().x(), rct2.top(),
                         rct2.center().x(), rct2.bottom())
        painter.setPen(pen)

        # draw two ellipses
        if True:
            diam = 5
            ymax = 5
            try:
                ymax = self.axisScaleDiv(Qwt.QwtPlot.yLeft).upperBound()
            except AttributeError:
                pass

            rect3 = QtCore.QRectF(rct)
            rect3.setLeft(xMap.transform(self.l_width) -  diam/2.0)
            rect3.setTop(yMap.transform(ymax/2) )
            rect3.setHeight(diam)
            rect3.setWidth(diam)
            painter.drawEllipse(rect3.toRect())

            rect3.setLeft(xMap.transform(self.r_width) -  diam/2.0)
            rect3.setWidth(diam)
            painter.drawEllipse(rect3.toRect())

    def update_all_curves(self, data, labels, ranges, mscore, intensity, assay_rt, show_legend):

        assert len(data) == len(labels)
        # This takes about 70% of the time
        self.create_curves(labels, ranges, show_legend)

        self.labels = []
        self.has_mscore = False
        if mscore is not None and len(ranges) < 2:
            self.labels = []
            self.has_mscore = True

            self.l_width = ranges[0][0]
            self.r_width = ranges[0][1]

            # create and add labels -> see drawCanvas
            try:
                mscore_txt = QtGui.QStaticText ("m_score %0.4g" % mscore);
                intensity_txt = QtGui.QStaticText ("Intensity %0.4g" % intensity);
                width_txt = QtGui.QStaticText ("PeakWidth %0.3fs" % (self.r_width - self.l_width) )
                assay_txt = QtGui.QStaticText ("Assay RT %0.3fs" % assay_rt)
                self.labels.extend([mscore_txt, intensity_txt, width_txt, assay_txt])
            except AttributeError:
                # old qt version?
                pass

        for d, curve in zip(data, self.curves):
            curve.setData( d[0], d[1] )

        self.replot()
        # self.get_plot().do_autoscale()
        self.clearZoomStack()

    def set_x_limits(self, xmin, xmax):
        self.setAxisScale(Qwt.QwtPlot.xBottom, xmin, xmax)

    def set_y_limits(self, ymin, ymax):
        self.setAxisScale(Qwt.QwtPlot.yLeft, ymin, ymax)

    def set_y_limits_auto(self, xmin, xmax):
        """
        Automatically set y limits based on the data in the plot and the
        visible range xmin, xmax.
        """
        allmin = []
        allmax = []
        for curve in self.curves:
            data = curve.data() # QwtArrayData
            filtered_data = [y for x,y in zip(data.xData(), data.yData()) if x>xmin and x<xmax]
            if len(filtered_data) == 0:
                continue
            allmin.append( min(filtered_data) )
            allmax.append( max(filtered_data) )
        if len(allmin) == 0 or len(allmax) == 0 : return
        self.set_y_limits(min(allmin), max(allmax) )

    def mousePressEvent(self, event):
        pass

    def mouseReleaseEvent(self, event):

        if 'lowerBound' in dir(self.axisScaleDiv(Qwt.QwtPlot.xBottom)):
            xmin = self.axisScaleDiv(Qwt.QwtPlot.xBottom).lowerBound()
            xmax = self.axisScaleDiv(Qwt.QwtPlot.xBottom).upperBound()
            ymin = self.axisScaleDiv(Qwt.QwtPlot.yLeft).lowerBound()
            ymax = self.axisScaleDiv(Qwt.QwtPlot.yLeft).upperBound()
        else:
            xmin = self.axisScaleDiv(Qwt.QwtPlot.xBottom).lBound()
            xmax = self.axisScaleDiv(Qwt.QwtPlot.xBottom).hBound()
            ymin = self.axisScaleDiv(Qwt.QwtPlot.yLeft).lBound()
            ymax = self.axisScaleDiv(Qwt.QwtPlot.yLeft).hBound()

        self.zoomChanged.emit(xmin, xmax, ymin, ymax)

