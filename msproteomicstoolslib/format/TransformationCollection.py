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
import msproteomicstoolslib.math.Smoothing as smoothing
import numpy

class LightTransformationData:
    """
    A lightweight data structure to store a transformation between retention times of multiple runs.
    """

    def __init__(self, ref=None):
        self.data = {} 
        self.trafo = {} 
        self.stdevs = {} 
        self.reference = ref

    def addTrafo(self, run1, run2, trafo, stdev=None):
      """
      Add transformation between two runs
      """
      d = self.trafo.get(run1, {})
      d[run2] = trafo
      self.trafo[run1] = d

      d = self.stdevs.get(run1, {})
      d[run2] = stdev
      self.stdevs[run1] = d

    def addData(self, run1, data1, run2, data2, doSort=True):
      """
      Add raw data for the transformation between two runs
      """
      # Add data from run1 -> run2 and also run2 -> run1
      assert len(data1) == len(data2)
      self._doAddData(run1, data1, run2, data2, doSort)
      self._doAddData(run2, data2, run1, data1, doSort)

    def _doAddData(self, run1, data1, run2, data2, doSort):
      if doSort and len(data1) > 0:
          data1, data2 = zip(*sorted(zip(data1, data2)))
      data1 = numpy.array(data1)
      data2 = numpy.array(data2)
      d = self.data.get(run1, {})
      d[run2] = (data1,data2)
      self.data[run1] = d

    #
    ## Getters
    #
    def getData(self, run1, run2):
        return self.data[run1][run2]

    def getTrafo(self, run1, run2):
        return self.trafo[run1][run2]

    def getStdev(self, run1, run2):
        return self.stdevs[run1][run2]

    def getTransformation(self, run1, run2):
        return self.trafo[run1][run2]

    def getReferenceRunID(self):
        return self.reference

class TransformationCollection():
    """A class to store a transformation between retention times of multiple runs.
    
    It allows to add transformation data (e.g. a pair of arrays which map
    coordinates from one RT space to the other). Once all data is added, one
    can initialize from the data:

    Compute a new transformation and write to file:
        # data1 = reference data (master) with ref_id
        # data2 = data to be aligned (slave) with current_id

    >>> tcoll = TransformationCollection()
    >>> tcoll.setReferenceRunID( ref_id )
    >>> tcoll.addTransformationData([data2, data1], current_id, ref_id )
    >>> tcoll.writeTransformationData( "outfile", current_id, ref_id)

    Read a set of transformations from files:

    >>> tcoll = TransformationCollection()
    >>> for filename in ["file1.tr", "file2.tr"]:
    >>>   tcoll.readTransformationData(filename)
    >>> tcoll.initialize_from_data(reverse=True)

    Compute a transformation:

    >>> norm_value = tcoll.getTransformation(orig_runid, ref_id).predict( [ value ] )[0]
    """

    def __init__(self):
      self.transformations = {}
      self.transformation_data = {}
      self._transformed_data = {}
      self._reference_run_id = None

    def getReferenceRunID(self):
        return self._reference_run_id

    def setReferenceRunID(self, value):
        self._reference_run_id = value

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

    def initialize_from_data(self, reverse=False, smoother="lowess", force=False):

        # use the data in self.transformation_data to create the trafos
        for s_from, darr in self.transformation_data.items():
            self.transformations[s_from] = {}
            import time
            for s_to, data in darr.items():
                start = time.time()
                if not self.getTransformedData(s_from, s_to) is None:
                    sm = smoothing.SmoothingInterpolation()
                    sm.initialize(data[0], self.getTransformedData(s_from, s_to))
                    self._addTransformation(sm, s_from, s_to)
                    if reverse: 
                        sm_rev = smoothing.SmoothingInterpolation()
                        sm_rev.initialize(self.getTransformedData(s_from, s_to), data[0])
                        self._addTransformation(sm_rev, s_to, s_from)
                else:
                    sm = smoothing.getSmoothingObj(smoother)
                    sm.initialize(data[0], data[1])
                    self.transformations[s_from][s_to] = sm
                    if reverse: 
                        sm_rev = smoothing.getSmoothingObj(smoother)
                        sm_rev.initialize(data[1], data[0])
                        self._addTransformation(sm_rev, s_to, s_from)
                print("Took %0.4fs to align %s against %s" % (time.time() - start, s_to, s_from))

    def addTransformationData(self, data, s_from, s_to):
      """ Add raw data points to the collection 

      Args:
          data(list(data_slave, data_master)) : two data two data vectors
              containing the raw data points from two runs. The first data
              vector is the master (reference) data and the second one is the
              slave (to be aligned).
          s_from(String): run ID of the slave (to be aligned) run
          s_to(String): run ID of the master (reference) run
      """
      assert isinstance(data, list)
      assert len(data) == 2
      assert isinstance(data[0], list) or isinstance(data[0], numpy.ndarray)
      assert isinstance(data[1], list) or isinstance(data[1], numpy.ndarray)
      d = self.transformation_data.get(s_from, {})
      d[s_to] = data
      self.transformation_data[s_from] = d

    def addTransformedData(self, data, s_from, s_to):
      """ Add transformed data points to the collection 

      The idea is to add the anchor points of s_from in the space of s_to so
      that one could compute the transformation using a simple linear transform.

      Args:
          data
          s_from(String): run ID of the slave (to be aligned) run
          s_to(String): run ID of the master (reference) run
      """
      assert isinstance(data, list)
      assert not self.getTransformationData(s_from, s_to) is None
      assert len(data) == len(self.getTransformationData(s_from, s_to)[0])
      d = self._transformed_data.get(s_from, {})
      d[s_to] = data
      self._transformed_data[s_from] = d

    def getTransformedData(self, s_from, s_to):
      try:
          return self._transformed_data[s_from][s_to]
      except KeyError:
          return None

    def getTransformationData(self, s_from, s_to):
      try:
          return self.transformation_data[s_from][s_to]
      except KeyError:
          return None

    def printTransformationData(self, s_from, s_to):
      r = self.getTransformationData(s_from, s_to)
      if r is None: return
      print("This data is able to transform from %s to %s" % (s_from, s_to))

    def writeTransformationData(self, filename, s_from, s_to):
      """Write the transformation s_from to s_to to a file.

      The header is either:
          #Transformation Null
          #Transformation Data "from_id" to "to_id" reference_id "ref_id"
      """
      r = self.getTransformationData(s_from, s_to)
      if r is None: 
          f = open(filename, "w")
          f.write("#Transformation Null\t%s" % s_from)
          f.close()
          return

      f = open(filename, "w")
      f.write("#Transformation Data\t%s\tto\t%s\treference_id\t%s\n" % (s_from, s_to, self._reference_run_id) )
      if s_from in self._transformed_data \
         and s_to in self._transformed_data[s_from]:
          tr = self._transformed_data[s_from][s_to]
          for a,b,c in zip(r[0],r[1], tr):
              f.write("%s\t%s\t%s\n" % (a,b,c) )
      else:
          for a,b in zip(r[0],r[1]):
              f.write("%s\t%s\n" % (a,b) )
      f.close()

    def readTransformationData(self, filename):
      """Read the transformation present in the file.

      The header is either:
          #Transformation Null
          #Transformation Data "from_id" to "to_id" reference_id "ref_id"
      """
      f = open(filename, "r")
      header = next(f).split("\t")
      if header[0].startswith( "#Transformation Null" ):
          # read the (or a) null transformation
          return
      s_from = header[1]
      s_to = header[3]
      if self._reference_run_id is None:
        self._reference_run_id = header[5].strip()
      assert self._reference_run_id == header[5].strip()
      data1 = []
      data2 = []
      transformed_data = []
      for line in f:
          d = line.split("\t")
          if len(d) != 2 and len(d) != 3: raise Exception("Cannot parse line '%s' in file %s" % (line, filename))
          data1.append(float(d[0]))
          data2.append(float(d[1]))
          if len(d) == 3:
            transformed_data.append(float(d[2]))
      # print("read data from %s to %s " %(s_from, s_to), [data1, data2])
      self.addTransformationData([data1, data2], s_from, s_to)
      if len(transformed_data) > 0: self.addTransformedData(transformed_data, s_from, s_to)

