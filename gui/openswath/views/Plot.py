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


from guiqwt.curve import CurveItem
from guiqwt.plot import CurvePlot, CurveDialog

from PyQt4 import QtGui 
from PyQt4.QtCore import Qt
from PyQt4 import QtCore

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
class GuiQwtMultiLinePlot(CurveDialog):
    """For the Curve window we could use a CurveDialog or a CurvePlot. 

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

    def create_curves(self, labels, this_range):

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

    def replot(self):
        self.get_plot().replot()

    def mouseReleaseEvent(self, event):
        xaxis_limits = self.get_plot().get_axis_limits("bottom")
        yaxis_limits = self.get_plot().get_axis_limits("left")
        self.zoomChanged.emit(xaxis_limits[0], xaxis_limits[1], yaxis_limits[0], yaxis_limits[1])

