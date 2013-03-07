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
    A small class to abstract some commonly used LaTeX commands.
"""

def prepare_float(myfloat, maxdigits, maxafter_comma=2):
    intpart = int( myfloat )
    floatpart = abs(myfloat - intpart)
    intstr = prepare_int(intpart, maxdigits)
    lenfloatpart = len(str(floatpart)) - maxafter_comma
    removedigits = -(lenfloatpart - 2)
    if removedigits > -1: return intstr + '.' + str(floatpart)[2:]
    return intstr + '.' + str(floatpart)[2:removedigits]

    #return intstr + '.' + str(floatpart)[2:-(lenfloatpart - 2)]

def prepare_int(myint, maxdigits):
    mystr = thousand_separate( myint)
    digits = maxdigits  - len(str(myint))
    dummy_nr = int( '1' * digits + str(myint) )
    padding = thousand_separate(dummy_nr)[:-len(mystr)]
    pad = '\\phantom{' + padding + '}'
    if digits == 0: return mystr
    return pad + mystr

def thousand_separate(myint):
    #we have to reverse to find every third digit
    tmp_mystr = ''
    for i,d in enumerate(reversed(str(myint))):
        if i % 3 == 0 and i != 0: tmp_mystr += ',\\'
        tmp_mystr += d
    #now go back to original
    ##mystr = ''
    ##for d in reversed(tmp_mystr):
    ##    mystr += d
    ##return mystr
    return reduce(lambda x,y: x+y, reversed(tmp_mystr))

def get_latex_row(arr):
    mystr = ''
    for item in arr:
        mystr += str( item )
        mystr += ' & '
    return mystr[:-2] + '\\\\'


latex_headers = """
\documentclass[pdftex, a4paper, landscape, 12pt, halfparskip, idxtotoc]{scrartcl}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}
\usepackage{array}
\usepackage{graphicx}
\usepackage{amsmath}
\DeclareGraphicsExtensions{.pdf,.jpg,.png}
\usepackage{url}
"""

def create_document(content, fname = 'result'):
    mydocument = latex_headers + r"\begin{document}"
    mydocument += content
    mydocument += "\end{document} "
    open( fname + '.tex', 'w').write(mydocument)
    os.system( 'pdflatex %s.tex' % fname )

def get_longtable(content, colnr, header, bold_caption, long_caption, label, cols=None):
    if cols is None: cols = "|" + "l| " * colnr
    headers = header.split( '&' )
    myheader = r"\multicolumn{1}{|c|}{\textbf{%s}}" % headers[0]
    for hea in headers[1:]:
        myheader += '\n' +  r"& \multicolumn{1}{c|}{\textbf{%s}}" % hea
    mytable = r"""
    \begin{longtable}{%(cols)s}
    \caption[%(c1)s]{\textbf{%(c1)s} %(c2)s}
    \label{%(label)s} \\
    \hline 
    %(head)s \\
    \hline 
    \endfirsthead
    \multicolumn{%(colnr)s}{c}
    {{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
    \hline 
    %(head)s \\
    \hline 
    \endhead
    \hline \multicolumn{%(colnr)s}{|r|}{{Continued on next page}} \\ \hline
    \endfoot
    \hline \hline
    \endlastfoot
    %(content)s
    \end{longtable}
    """ % { 'c1' : bold_caption, 'c2' : long_caption, 'label' : label,
           'colnr' : colnr, 'head' : myheader, 'content' : content, 'cols' : cols }
    return mytable

def get_table(content, colnr, header, bold_caption, long_caption, label, cols=None):
    if cols is None: cols = "c " * colnr
    mytable = r"""
    \begin{table}[h]
    \centering
    \caption[%(c1)s]{\textbf{%(c1)s} %(c2)s}
    \label{tab:%(label)s}
    \begin{tabular}{ %(cols)s }
    %(tabletop)s \\
    %(content)s
    \end{tabular}
    \end{table}
    """ % { 'c1' : bold_caption, 'c2' : long_caption, 'label' : label,
           'cols' : cols, 'tabletop' : header, 'content' : content }
    return mytable



""
#testcases
def test():
    assert get_latex_row( [5, 'a', 99]) == '5 & a & 99 \\\\'

    testing = [200, 5000, 60000, 800000, 9000000]
    testing2 = [20.10, 5.034, -160.400, 0.856780, 9.000000]
    result = ''
    for myc in range(5):
        number  = prepare_int(testing[myc], 7)
        number2  = prepare_float(testing2[myc], 3)
        result += get_latex_row( [myc, number, number2 ] ) + '\n'

    #TODO handle the minus correctly
    assert result == '0 & \\phantom{1\\,111\\,}200 & \\phantom{1}20. \\\\\n1 & \\phantom{1\\,11}5\\,000 & \\phantom{11}5.03 \\\\\n2 & \\phantom{1\\,1}60\\,000 & \\phantom{}-\\,160. \\\\\n3 & \\phantom{1\\,}800\\,000 & \\phantom{11}0.85 \\\\\n4 & 9\\,000\\,000 & \\phantom{11}9. \\\\\n'

    assert prepare_float( 35.82 , 2, 2) == '35.82'
    assert prepare_float( 35.82 , 4, 2) == '\\phantom{1\\,1}35.82'
    assert prepare_float( 35.82 , 2, 3) == '35.82'

