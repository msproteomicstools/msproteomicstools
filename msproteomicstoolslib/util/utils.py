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
from __future__ import print_function

"""
Doc :
    A random collection of useful and not so useful functions and objects.
"""

import sys

class stringp:
        def contains(theString, theQueryValue):
            return theString.find(theQueryValue) > -1
        
        def delSequenceDots(seq):
            if seq.find('.')==2:
                sequenceNoDots = seq[3:-3]
            else:
                sequenceNoDots = seq

            return sequenceNoDots

class Math:

    def median(y):
        z = len(y)
        if not z%2: return (y[(z/2)-1] + y[z/2]) / 2
        else: return y[z/2]        


class Lists:
    
    def getNotFound(self, bigList, smallList):
        """Returns the elements of smallList not contained in bigList"""
        
        if len(bigList)<len(smallList):
            print("bad lists provided!")
            sys.exit()
            
        resultingList = list(smallList)

        for elem in smallList:
            if elem in bigList:
                resultingList.remove(elem)
        
        return resultingList
    
    def getFound(self, list1, list2):
        """Returns the common elements of two lists"""

        bigList= []
        smallList = []
        commonList = []
        
        if len(list1)>=len(list2):
            bigList = list1
            smallList = list2
        else:
            bigList = list2
            smallList = list1
        
        for elem in smallList:
            if elem in bigList:
                commonList.append(elem)
        
        return commonList
        
        
    
    def readPeptideList(self, file):
        import csv
        mylist = []
        csvfile = open(file,"rb")
        reader = csv.reader(csvfile, delimiter='\t',quotechar=' ')
        
        for row in reader:
            for k in row:
                mylist.append(k.strip())
                        
        csvfile.close()
        
        return mylist
  
    


def attach_method(fxn, instance, object):
    """We add the function fxn to the instance of a certain object.

    The function can access all attributes stored as self; the first argument
    of the function should be self like in all normal functions of a class.
    This is helpful if you want to add a function after creation. example

    def fxn(self, ...):
        self.attribute
    """
    import types
    f = types.MethodType(fxn, instance, object)
    exec('instance.' + fxn.__name__ + ' = f')

def unique(seq): 
    # order preserving
    checked = []
    for e in seq:
        if e not in checked:
            checked.append(e)
    return checked

class db_table:

    def __init__(self, cursor):
        self.c = cursor

    def read(self, query, arr = []):
        self.result = self.c.execute(query, arr)
        #get the coloumn description out of it as hash
        #element zero is the name of the coloumn
        self.desc = {}
        for d in range(len(self.c.description)):
            self.desc[ self.c.description[d][0] ] = d

    def fetchone(self):
        self.lastRow = self.c.fetchone()
        return self.lastRow

    def get_col(self, name):
        return self.lastRow[ self.desc[ name ] ]

    def row(self, row, name):
        return row[ self.desc[ name ] ]

    def fetchall_groupBy(self, name):
        resHash = {}
        while True:
            row = self.c.fetchone()
            if row == None: break
            gr_name = row[ self.desc[ name ] ]
            if gr_name in resHash: resHash[gr_name].append(row)
            else: resHash[gr_name] = [row]
        return resHash

    def rows(self):
        while True:
            row = self.c.fetchone()
            if row == None: break 
            yield row

    def rows_groupby(self, name):
        """Iterator -- attention this assumes that the rows are sorted!"""
        start = True
        current = None
        stack = []
        while True:
            row = self.fetchone()
            if row == None: break
            gr_name = row[ self.desc[ name ] ]
            if gr_name == current: stack.append( row )
            else: 
                if len(stack) > 0: yield stack
                stack = [ row ]
                current = gr_name

class csv_table:

    def __init__(self, filename, delimiter = '\t', header = True):
        self.filename = filename
        self.delimiter = delimiter
        self.header_true = header
        self.header = {}

    def read(self):
        import csv
        r = csv.reader(open(self.filename), delimiter=self.delimiter)
        self.rows = []
        if self.header_true: 
            self.header_arr = r.next()
            for i, head in enumerate(self.header_arr):
                self.header[ head ] = i
        for row in r:
            self.rows.append(row)

    def rows(self):
        import csv
        r = csv.reader(open(self.filename), delimiter=self.delimiter)
        if self.header_true: 
            self.header_arr = r.next()
            for i, head in enumerate(self.header_arr):
                self.header[ head ] = i
        for row in r:
            yield row

    def get_coloumn(self, name):
        i = self.header[ name ]
        if name in self.header:
            res = [ tmp[i] for tmp in self.rows]
        return res

