#!/usr/bin/python
# -*- coding: utf-8  -*-
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

import os


"""
Doc :
    A small class to abstract some commonly used Gnuplot commands.
"""

class Gnuplot:
    def __init__(self, nocolor=False):
        self.header = "set terminal postscript color enhanced solid\n"
        if nocolor: self.header = "set terminal postscript enhanced\n"
        #thin border
        self.header += "set border 31 linewidth .3\n"
        self.body = ''
        self.has_title = False
        self.output = 'set output "test.eps"'
    def set_x_y_label(self, xlabel, ylabel):
        self.body += 'set xlabel "%s"\n' % xlabel
        self.body += 'set ylabel "%s"\n' % ylabel
    def set_output(self, out):
        self.output = 'set output "%s"\n' % out
    def set_input(self, input):
        self.input = input
    def set_nokey(self):
        self.body += 'set nokey\n' 
    def add_to_body(self, addition):
        self.body += addition + '\n'
    def set_title(self, title):
        self.body += 'set title "%s"\n' % title
        self.has_title = True
    def draw_points(self, lw = 2, keep_data = False):
        if self.has_title: notitle = ''
        else: notitle = 'notitle'
        self.plot = 'plot "%s" %s lt -1 lw %s\n' % (self.input, notitle, lw)
        tmp_filename = 'gnuplot.plt'
        #
        self.gnuplot = self.header + self.body + self.output + self.plot
        f = open( tmp_filename, 'w')
        f.write( self.gnuplot )
        f.close()
        os.system( "gnuplot %s" % tmp_filename )
        if not keep_data: os.system( "rm %s" % tmp_filename )
    def draw_boxes(self, lw = 2, keep_data = False):
        if self.has_title: notitle = ''
        else: notitle = 'notitle'
        self.plot = 'plot "%s" %s with boxes lt -1 lw %s\n' % (self.input, notitle, lw)
        tmp_filename = 'gnuplot.plt'
        #
        self.gnuplot = self.header + self.body + self.output + self.plot
        f = open( tmp_filename, 'w')
        f.write( self.gnuplot )
        f.close()
        os.system( "gnuplot %s" % tmp_filename )
        if not keep_data: os.system( "rm %s" % tmp_filename )
    def draw(self, cmd='with points lt -1 lw 2', keep_data = False, plot=None):
        if self.has_title: notitle = ''
        else: notitle = 'notitle'
        self.plot = 'plot "%s" %s %s\n' % (self.input, notitle, cmd)
        if not plot is None: self.plot = plot
        tmp_filename = '/tmp/gnuplot.plt'
        #
        self.gnuplot = self.header + self.body + self.output + self.plot
        f = open( tmp_filename, 'w')
        f.write( self.gnuplot )
        f.close()
        os.system( "gnuplot %s" % tmp_filename )
        if not keep_data: os.system( "rm %s" % tmp_filename )
    @staticmethod
    def draw_boxes_from_data(data, output = 'test.eps', 
            xlabel = '', ylabel = '', tmp_csv = 'tmp.csv', title=None, 
            keep_data = False):
        #
        #
        n = data[0]; h = data[1]
        f = open( tmp_csv , 'w')
        for nn, hh in zip(n,h):
            f.write( '%s %s\n' % (hh, nn)  )
        f.close()
        gnu = Gnuplot()
        gnu.set_x_y_label( xlabel, ylabel )
        gnu.set_output( output )
        gnu.set_input(  tmp_csv )
        if not title == None: gnu.set_title( title ); gnu.set_nokey()
        gnu.draw_boxes(keep_data=keep_data)
        if not keep_data: os.system( 'rm %s' % tmp_csv)
        return gnu
    @staticmethod
    def draw_points_from_data(data, output = 'test.eps', 
            xlabel = '', ylabel = '', tmp_csv = 'tmp.csv', title=None, 
            keep_data = False):
        #
        #
        n = data[0]; h = data[1]
        f = open( tmp_csv , 'w')
        for nn, hh in zip(n,h):
            f.write( '%s %s\n' % (hh, nn)  )
        f.close()
        gnu = Gnuplot()
        gnu.set_x_y_label( xlabel, ylabel )
        gnu.set_output( output )
        gnu.set_input(  tmp_csv )
        if not title == None: gnu.set_title( title ); gnu.set_nokey()
        gnu.draw_points()
        if not keep_data: os.system( 'rm %s' % tmp_csv)
        return gnu
    @staticmethod
    def draw_from_data(data, cmd='with points lt -1 lw 2', output = 'test.eps', 
            xlabel = '', ylabel = '', tmp_csv = 'tmp.csv', title=None, 
            keep_data = False):
        #
        #
        n = data[0]; h = data[1]
        f = open( tmp_csv , 'w')
        for nn, hh in zip(n,h):
            f.write( '%s %s\n' % (hh, nn)  )
        f.close()
        gnu = Gnuplot()
        gnu.set_x_y_label( xlabel, ylabel )
        gnu.set_output( output )
        gnu.set_input(  tmp_csv )
        if not title == None: gnu.set_title( title ); gnu.set_nokey()
        gnu.draw(cmd)
        if not keep_data: os.system( 'rm %s' % tmp_csv)
        return gnu

    @staticmethod
    def draw_from_multiple_data(data, cmd='with points lt -1 lw 2', output = 'test.eps', 
            xlabel = '', ylabel = '', tmp_csv = '/tmp/htmp', title=None, datatitles=None,
            keep_data = False, body_add='', nocolor=False, addls=False):
        #
        #
        for i,(n,h) in enumerate(data): 
            filename = tmp_csv + '%s.csv' % i
            f = open( filename , 'w')
            for nn, hh in zip(n,h):
                f.write( '%s %s\n' % (hh, nn)  )
            f.close()

        gnu = Gnuplot(nocolor)
        gnu.set_x_y_label( xlabel, ylabel )
        gnu.set_output( output )
        gnu.add_to_body(body_add)
        gnu.plot = 'plot '
        notitle = ''
        for i,(n,h) in enumerate(data): 
            filename = tmp_csv + '%s.csv' % i
            if datatitles: notitle = 'title "%s"' %  datatitles[i]
            if addls: gnu.plot += ' "%s" %s %s ls %s, \\\n' % (filename, notitle, cmd, i)
            else: gnu.plot += ' "%s" %s %s, \\\n' % (filename, notitle, cmd)

        gnu.plot = gnu.plot[:-4]
        if not title == None: gnu.set_title( title )
        gnu.input = ''
        gnu.draw(cmd, plot=gnu.plot, keep_data=keep_data)
        if not keep_data: 
            for i,(n,h) in enumerate(data): 
                filename = tmp_csv + '%s.csv' % i
                os.system( 'rm %s' % filename)
        return gnu

