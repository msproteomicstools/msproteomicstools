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

import numpy

def get_smooting_operator():
  try:
    import rpy2.robjects as robjects
    return SmoothingR()
  except ImportError:
    pass
  # Now try the Python
  try:
    from scikits import datasmooth as ds
    scikits_present = True
    return SmoothingPy()
  except ImportError:
    pass
  return None

class SmoothingR:
    """Class to smooth data using the smooth.spline function from R

     data1 = c(5,7,8,9,10,15,7.1,6)
     data2 = c(4,7,9,11,11,14,7.1,6.5)
     data1 = sort(data1)
     data2 = sort(data2)
     smooth.model = smooth.spline(data1,data2,cv=T)
     data2_pred = predict(smooth.model,data2)$y 
     [1]  2.342662  6.615797  7.292613  7.441842 10.489440 11.858406 11.858406
     [8] 13.482255
     plot(data1, data2)
     lines(data1, data2_pred, col="blue")


    # In python ...
    import rpy2.robjects as robjects
    # uses python-rpy2
    data1 = [5,7,8,9,10,15,7.1,6]
    data2 = [4,7,9,11,11,14,7.1,6.5]
    rdata1 = robjects.FloatVector(data1)
    rdata2 = robjects.FloatVector(data2)
    spline = robjects.r["smooth.spline"]
    sm = spline(data1,data2,cv=T)
    predict = robjects.r["predict"]
    predicted_data = predict(sm, rdata2)
    numpy.array(predicted_data[1])
    array([  2.34266247,   7.2926131 ,  10.48943975,  11.85840597,
            11.85840597,  13.48225519,   7.44184246,   6.61579704])

    """

    def __init__(self):
        try:
          import rpy2.robjects as robjects
        except ImportError:
            print "==================================="
            print "rpy2 package, please install it first\n (see https://pypi.python.org/pypi/rpy2/)." 
            print "==================================="

    def initialize(self, data1, data2):
        import rpy2.robjects as robjects
        rdata1 = robjects.FloatVector(data1)
        rdata2 = robjects.FloatVector(data2)
        spline = robjects.r["smooth.spline"]
        self.sm = spline(data1,data2,cv=True)

    def predict(self, xhat):
        import rpy2.robjects as robjects
        rxhat = robjects.FloatVector(xhat)
        predict = robjects.r["predict"]
        predicted_data = predict(self.sm, rxhat)
        predicted_result = numpy.array(predicted_data[1]).tolist()
        return predicted_result

class SmoothingPy:
    """Smoothing of 2D data using generalized crossvalidation

    Will call _smooth_spline_scikit internally but only at a few select
    points. It then uses the generated smoothed spline to construct an
    interpolated spline on which then the xhat data is evaluated.
    """

    def __init__(self):
        import operator
        try:
          from scikits import datasmooth as ds
        except ImportError:
            print "==================================="
            print "Cannot import the module datasmooth from scikits, \nplease download it from https://github.comtickel/scikit-datasmooth.git"
            print "==================================="

    def de_duplicate_array(self, arr):
        arr_fixed = [] 
        duplications = []
        i = 0
        while i < numpy.size(arr)-1:
            dupl = 1
            if arr[i] == arr[i+1]: 
                while i < numpy.size(arr)-1 and arr[i] == arr[i+1] : 
                    i += 1
                    dupl += 1 # indices.append(i)
            arr_fixed.append(arr[i])
            duplications.append(dupl)
            i += 1
        return arr_fixed, duplications

    def re_duplicate_array(self, arr_fixed, duplications):
        # arr = [0, 0, 5, 6, 6, 7, 8, 8]
        result = [] 
        for val, dupl in zip(arr_fixed, duplications):
            result.extend( [val for i in range(dupl) ] )
        return result

    def _smooth_spline_scikit(self, data1, data2, xhat=None, fixNonMonotonous=False):
        """Smoothing of 2D data using generalized crossvalidation

        Will return the evaluated data at the points xhat (or if xhat is empty,
        at the points of data1). Note that the result _will_ depend on xhat
        since the optimization function will try to maximize the smoothness of
        the line generated by (xhat,yhat).

        Do not call this function with a large set of (or too densily spaced)
        evaluation data. Rather use the wrap function.

        uses datasmooth from https://github.comtickel/scikit-datasmooth.git
        """
        try:
          from scikits import datasmooth as ds
        except ImportError:
            print "==================================="
            print "Cannot import the module datasmooth from scikits, \nplease download it from https://github.comtickel/scikit-datasmooth.git"
            print "==================================="
            import sys; sys.exit(1)
        import operator

        x = numpy.array(data1)
        y = numpy.array(data2)
        if xhat is None:
            xhat = numpy.array(x)
        else:
            xhat = numpy.array(xhat)

        # Step 1: get the indices of the original xhat
        tmp_xhat = sorted(enumerate(xhat), key=operator.itemgetter(1))
        xhat_sorted = numpy.array([t[1] for t in tmp_xhat])
        xhat_indices = [t[0] for t in tmp_xhat]

        # fix if not monotonous increasing... 
        if fixNonMonotonous:
            xhat_sorted,duplications = self.de_duplicate_array(xhat_sorted)
            xhat_sorted = numpy.array(xhat_sorted)

        # Step 2: Execute the call to smooth the data
        # throws memory error for large data -> use the wrapper 
        yhat_sorted,lmbd = ds.smooth_data(x,y,xhat=xhat_sorted)

        # Step 3.1 re-insert duplicated values
        if fixNonMonotonous:
            yhat_sorted = self.re_duplicate_array(yhat_sorted, duplications)

        # Step 3: re-order the values using the original indices from before
        yhat = [None for i in range(numpy.size(xhat))]
        for i in range(numpy.size(xhat)):
            yhat[ xhat_indices[i] ] = yhat_sorted[i]

        return yhat

        # plot(x,y,'ow',xhat_sorted,yhat_sorted,'-b', x,yhat, "or")

    def initialize(self, data1, data2, Nhat=200, xmin=None, xmax=None):
        from scipy.interpolate import InterpolatedUnivariateSpline

        # Create a spline with Nhat points
        if xmin is None: xmin = numpy.min(data1)
        if xmax is None: xmax = numpy.max(data1)
        xh = numpy.linspace(xmin-0.1,xmax+1.1,Nhat)
        yhat_small = self._smooth_spline_scikit(data1, data2, xh)

        # Now use that (small) smoothed spline to create a univariate spline 
        self.ius = InterpolatedUnivariateSpline(xh, yhat_small)

    def predict(self, xhat):
        xhat = numpy.array(xhat)
        yhat_new = self.ius(xhat)
        return yhat_new

    # TODO remove
    def _smooth_scikit_legacy(self, data1, data2, xhat, Nhat=200):
        xhat = numpy.array(xhat)
        xmin = numpy.min(xhat)
        xmax = numpy.max(xhat)

        s = SmoothingPy()
        s.initialize(data1, data2, Nhat, xmin = xmin, xmax = xmax)
        yhat_new = s.predict(xhat)
        return yhat_new
