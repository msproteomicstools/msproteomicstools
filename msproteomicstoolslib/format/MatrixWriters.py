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
$Authors: Hannes Roest, Lorenz Blum$
--------------------------------------------------------------------------
"""

from __future__ import print_function
import csv

try:
    import xlsxwriter
    import xlwt
except ImportError:
    # Python 3
    import xlwt


def getwriter(matrix_outfile):
    """
    Factory function to get the correct writer depending on the file ending

    Args:
        matrix_outfile(str): Filename of output - used to determine output format. Valid formats are .xlsx .xls .csv or .tsv
    """
    if matrix_outfile.endswith("xls"):
        matrix_writer = XlsWriter(matrix_outfile)
    elif matrix_outfile.endswith("xlsx"):
        matrix_writer = XlsxWriter(matrix_outfile)
    elif matrix_outfile.endswith("tsv"):
        matrix_writer = CsvWriter(matrix_outfile, delim="\t")
    elif matrix_outfile.endswith("csv"):
        matrix_writer = CsvWriter(matrix_outfile, delim=",")
    else:
        raise Exception("Unknown matrix extension, must be .xlsx .xls .csv or .tsv")

    return matrix_writer

class IWriter():
    """
    Interface. you need to implement init, write, newline and del
    """
    def __init__(self, outfile, delim=None):
        raise NotImplementedError

    def write(self, entry, color=None):
        raise NotImplementedError

    def newline(self):
        raise NotImplementedError

    def __del__(self):
        raise NotImplementedError


class CsvWriter(IWriter):
    def __init__(self, outfile, delim="\t"):
        self.outfile = outfile
        self.f = open(outfile, "w")
        self.tsv = csv.writer(self.f, delimiter=delim)
        self.line = []

    def write(self, entry, color="ignored"):
        self.line.append(entry)

    def newline(self):
        self.tsv.writerow(self.line)
        self.line = []

    def __del__(self):
        print("Closing ", self.outfile)
        self.f.close()


class XlsWriter(IWriter):
    def __init__(self, outfile, delim='ignored'):
        self.outfile = outfile
        self.workbook = xlwt.Workbook()
        self.worksheet = self.workbook.add_sheet('0')
        self.col = 0
        self.row = 0
        self.colors = {'r': xlwt.easyxf('font: color red;'), 'b': xlwt.easyxf('font: color blue;'), 'd': xlwt.easyxf()}

    def write(self, entry, color='d'):
        if color == 'd':
            self.worksheet.write(self.row, self.col, entry)
        else:
            thestyle = self.colors[color]
            self.worksheet.write(self.row, self.col, entry, thestyle)
        self.col += 1

    def newline(self):
        self.row += 1
        self.col = 0

    def __del__(self):
        print("Writing out", self.outfile)
        self.workbook.save(self.outfile)


class XlsxWriter(IWriter):
    def __init__(self, outfile, delim='ignored'):
        self.outfile = outfile
        self.workbook = xlsxwriter.Workbook(outfile)
        self.worksheet = self.workbook.add_worksheet()
        self.colors = {'r': 'red', 'b': 'blue', 'd': 'black'}
        self.col = 0
        self.row = 0

    def write(self, entry, color='d'):
        if 'nan' in str(entry):
            entry = str(entry)
        if color == 'd':
            self.worksheet.write(self.row, self.col, entry)
        else:
            fmt = self.workbook.add_format()
            fmt.set_font_color(self.colors[color])
            self.worksheet.write(self.row, self.col, entry, fmt)
        self.col += 1

    def newline(self):
        self.row += 1
        self.col = 0

    def __del__(self):
        print("Writing out", self.outfile)
        self.workbook.close()
