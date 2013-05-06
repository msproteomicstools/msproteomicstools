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

import msproteomicstoolslib.math.Smoothing as smoothing

class TransformationCollection():
    """A class to store information about multiple transformations between two
    runs.
    
    It allows to add transformation data (e.g. a pair of arrays which map
    coordinates from one RT space to the other). Once all data is added, one
    can initialize from the data:

    Example:

    >>> transformation_collection_ = TransformationCollection()
    >>> for filename in ["file1.tr", "file2.tr"]:
    >>>   transformation_collection_.readTransformationData(filename)
    >>> transformation_collection_.initialize_from_data(reverse=True)

    """
    def __init__(self):
      self.transformations = {}
      self.transformation_data = {}
      self.reference_run_id = None

    def _addTransformation(self, trafo, s_from, s_to):
      d = self.transformations.get(s_from, {})
      d[s_to] = trafo
      self.transformations[s_from] = d

    def getTransformation(self, s_from, s_to):

      if s_from == s_to: 
          # null smoothing
          return smoothing.SmoothingNull()

      try:
          return self.transformations[s_from][s_to]
      except KeyError:
          return None

    def initialize_from_data(self, reverse=False):
        # use the data in self.transformation_data to create the trafos
        for s_from, darr in self.transformation_data.iteritems():
            self.transformations[s_from] = {}
            for s_to, data in darr.iteritems():
                sm = smoothing.get_smooting_operator()
                sm.initialize(data[0], data[1])
                self.transformations[s_from][s_to] = sm

                if reverse: 

                    sm_rev = smoothing.get_smooting_operator()
                    sm_rev.initialize(data[1], data[0])
                    self._addTransformation(sm_rev, s_to, s_from)

    def addTransformationData(self, data, s_from, s_to):
      assert isinstance(data, list)
      assert len(data) == 2
      assert isinstance(data[0], list)
      assert isinstance(data[1], list)
      d = self.transformation_data.get(s_from, {})
      d[s_to] = data
      self.transformation_data[s_from] = d

    def getTransformationData(self, s_from, s_to):
      try:
          return self.transformation_data[s_from][s_to]
      except KeyError:
          return None

    def printTransformationData(self, s_from, s_to):
      r = self.getTransformationData(s_from, s_to)
      if r is None: return
      print "This data is able to transform from %s to %s" % (s_from, s_to)

    def writeTransformationData(self, filename, s_from, s_to):
      r = self.getTransformationData(s_from, s_to)
      if r is None: 
          f = open(filename, "w")
          f.write("#Transformation Null\t%s" % s_from)
          f.close()
          return

      f = open(filename, "w")
      f.write("#Transformation Data\t%s\tto\t%s\treference_id\t%s\n" % (s_from, s_to, self.reference_run_id) )
      for a,b in zip(r[0],r[1]):
          f.write("%s\t%s\n" % (a,b) )
      f.close()

    def readTransformationData(self, filename):
      f = open(filename, "r")
      header = f.next().split("\t")
      if header[0].startswith( "#Transformation Null" ):
          # read the (or a) null transformation
          return
      s_from = header[1]
      s_to = header[3]
      if self.reference_run_id is None:
        self.reference_run_id = header[5].strip()
      assert self.reference_run_id == header[5].strip()
      data1 = []
      data2 = []
      for line in f:
          d = line.split("\t")
          data1.append(float(d[0]))
          data2.append(float(d[1]))
      # print "read data from %s to %s " %(s_from, s_to), [data1, data2]
      self.addTransformationData([data1, data2], s_from, s_to)
