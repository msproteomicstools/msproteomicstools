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
import numpy
import os
import csv
import math
import random
from numpy import median, absolute

def get_smooting_operator(use_scikit=False, use_linear=False, use_external_r = False, tmpdir=None):
  if use_linear: 
      return SmoothingLinear()
  try:
    if use_scikit: 
        raise ImportError
    if use_external_r: 
        return SmoothingRExtern(tmpdir)
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
    print("No smoothing operator is available, please install either rpy2 or scikits with datasmooth.")
  return None

def getSmoothingObj(smoother, topN=5, max_rt_diff=30, min_rt_diff=0.1, removeOutliers=False, tmpdir=None):
    if smoother == "diRT":
        return SmoothingNull()
    elif smoother == "linear":
        return SmoothingLinear()
    elif smoother == "splineR":
        return SmoothingR()
    elif smoother == "splineR_external":
        return SmoothingRExtern()
    elif smoother == "splinePy":
        return SmoothingPy()
    elif smoother == "lowess":
        return LowessSmoothingStatsmodels()
    elif smoother == "lowess_biostats":
        return LowessSmoothingBiostats()
    elif smoother == "lowess_statsmodels":
        return LowessSmoothingStatsmodels()
    elif smoother == "lowess_cython":
        return LowessSmoothingCyLowess()
    elif smoother == "nonCVSpline":
        return UnivarSplineNoCV()
    elif smoother == "CVSpline":
        return UnivarSplineCV()
    elif smoother == "Earth":
        return SmoothingEarth()
    elif smoother == "WeightedNearestNeighbour":
        return WeightedNearestNeighbour(topN, max_rt_diff, min_rt_diff, removeOutliers)
    elif smoother == "SmoothLLDMedian":
        return SmoothLLDMedian(topN, max_rt_diff, min_rt_diff, removeOutliers)
    elif smoother == "None":
        return SmoothingNull()
    else:
        raise Exception("Unknown smoothing method: " + smoother)

class SmoothingR:
    """Class to smooth data using the smooth.spline function from R
    
    This is equivalent to the following code::

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


    Doing the same thing in Python ::

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
            print("===================================")
            print("rpy2 package, please install it first\n (see https://pypi.python.org/pypi/rpy2/)." )
            print("===================================")

    def initialize(self, data1, data2):
        import rpy2.robjects as robjects
        robjects.r['options'](warn=-1)
        rdata1 = robjects.FloatVector(data1)
        rdata2 = robjects.FloatVector(data2)
        spline = robjects.r["smooth.spline"]
        self.sm = spline(data1,data2,cv=True)

    def predict(self, xhat):
        import rpy2.robjects as robjects
        robjects.r['options'](warn=-1)
        rxhat = robjects.FloatVector(xhat)
        predict = robjects.r["predict"]
        predicted_data = predict(self.sm, rxhat)
        predicted_result = numpy.array(predicted_data[1]).tolist()
        return predicted_result

class SmoothingRExtern:
    """Class to smooth data using the smooth.spline function from R (extern system call)
    """

    def __init__(self, TMPDIR="/tmp/"):
        if TMPDIR is None:
            raise Exception("Tempdir needs to be set (cannot be none)")
        self.TMPDIR = TMPDIR

    def initialize(self, data1, data2):
        prediction_data = data2
        arr = self.predict_R_(data1, data2, prediction_data, self.TMPDIR)

        # Internally then use interpolation to actually predict
        self.internal_interpolation = SmoothingInterpolation()
        self.internal_interpolation.initialize(arr[:,0], arr[:,1])

    def predict_R_(self, data1, data2, predict_data, TMPDIR):
        fname = TMPDIR + "/datafile_feature_align_%s" % int(random.random() * 100000)
        fname_pred = TMPDIR + "/datafile_feature_align_%s" % int(random.random() * 100000)
        fname_out = TMPDIR + "/datafile_feature_align_%s" % int(random.random() * 100000)
        Rscript = TMPDIR + "/tmp"+str(int(random.random() * 100000))+".R"

        # Input file with datapoints
        fh = open(fname, "w")
        fh.write("data1\tdata2\n")
        for d1, d2 in zip(data1, data2):
            fh.write("%s\t%s\n" % (d1,d2))
        fh.close()

        # File which datapoints to predict
        fh = open(fname_pred, "w")
        fh.write("predict_on\n")
        for d2 in predict_data:
            fh.write("%s\n" % (d2))
        fh.close()

        fh = open(Rscript, "w")
        fh.write( """
        args <- commandArgs(trailingOnly = TRUE)
        # trailingOnly=TRUE means that only arguments after --args are returned
        infile <- args[1]
        predictionfile <- args[2]
        outfile <- args[3]
        df = read.csv(infile,header=TRUE, sep="\t")
        predict_on = read.csv(predictionfile,header=TRUE, sep="\t")
        sm = smooth.spline(df$data1,df$data2,cv=T)
        prediction = predict(sm,predict_on[,1])$y
        write.table(data.frame(predict_on, prediction), outfile, sep="\t", row.names=FALSE)
        """)
        fh.close()

        # Execute command
        cmd = "R --slave --args %s %s %s < %s" % (fname, fname_pred, fname_out, Rscript) 
        os.system(cmd)

        try:
            r = csv.reader(open(fname_out), delimiter="\t")
            next(r)
            arr = numpy.array([ (float(line[0]),float(line[1])) for line in r ])
        except IOError as e:
            print("Something went wrong, I cannot find the file at ", fname_out)
            print("Debug output:")
            print("Input data length (d1, d2):", len(data1), len(data2))
            print("Input data length:", len(predict_data))
            print("Temporary directory :", TMPDIR)
            raise e

        # Cleanup
        os.system("rm %s" % fname)
        os.system("rm %s" % fname_pred)
        os.system("rm %s" % fname_out)
        os.system("rm %s" % Rscript)

        return arr

    def predict(self, xhat):
        return self.internal_interpolation.predict(xhat)

class SmoothingNull:
    """Null smoother that performs a null operation """

    def __init__(self):
        pass

    def initialize(self, data1, data2):
        pass

    def predict(self, xhat):
        return xhat

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
            print("===================================")
            print("Cannot import the module datasmooth from scikits, \nplease download it from https://github.comtickel/scikit-datasmooth.git")
            print("===================================")

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
            print("===================================")
            print("Cannot import the module datasmooth from scikits, \nplease download it from https://github.comtickel/scikit-datasmooth.git")
            print("===================================")
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
        return list(yhat_new)

    # TODO remove
    def _smooth_scikit_legacy(self, data1, data2, xhat, Nhat=200):
        xhat = numpy.array(xhat)
        xmin = numpy.min(xhat)
        xmax = numpy.max(xhat)

        s = SmoothingPy()
        s.initialize(data1, data2, Nhat, xmin = xmin, xmax = xmax)
        yhat_new = s.predict(xhat)
        return yhat_new

class LowessSmoothingBase(object):
    """Smoothing using Lowess smoother and then interpolate on the result
    """

    def __init__(self):
        pass

    def initialize(self, data1, data2):
        result = self._initialize(data1, data2)

        self.internal_interpolation = SmoothingInterpolation()
        self.internal_interpolation.initialize(result[0], result[1])

    def predict(self, xhat):
        return self.internal_interpolation.predict(xhat)

class LowessSmoothingBiostats(LowessSmoothingBase):
    """Smoothing using Lowess smoother and then interpolate on the result
    """

    def __init__(self):
        pass

    def _initialize(self, data1, data2):
        try:
            from Bio.Statistics.lowess import lowess
        except ImportError:
            print("===================================")
            print("Cannot import the module lowess from Biopython, \nplease install 'biopython' from https://pypi.python.org/pypi/biopython")
            print("===================================")

        old_settings = numpy.seterr(all='ignore')

        result = lowess(numpy.array(data1), numpy.array(data2), f=0.1, iter=3)
        if all([math.isnan(it) for it in result]):
            # Try standard paramters
            result = lowess(numpy.array(data1), numpy.array(data2))

        numpy.seterr(**old_settings)
        return data1, result

class LowessSmoothingStatsmodels(LowessSmoothingBase):
    """Smoothing using Lowess smoother and then interpolate on the result

    statsmodels now also has fast Cython lowess, see https://github.com/statsmodels/statsmodels/pull/856

    This faster lowess should be in version 0.5.0 of statsmodels (anaconda
    currently has version 0.6.0). However, Ubuntu only has version 0.5.0 from
    14.04 onwards, so be careful.

	frac: float
        Between 0 and 1. The fraction of the data used
        when estimating each y-value.
        it: int
        The number of residual-based reweightings
        to perform.

    """

    def __init__(self):
        pass

    def _initialize(self, data1, data2):
        try:
            import statsmodels.api as sm
            lowess = sm.nonparametric.lowess
        except ImportError:
            print("===================================")
            print("Cannot import the module lowess from 'statsmodels', \nplease install the Python package 'statsmodels'")
            print("===================================")

        # NOTE: delta parameter is only available from statsmodels > 0.5.0
        delta = (max(data1) - min(data1)) * 0.01
        frac = 0.1
        
        if len(data1) < 100:
            frac = 1.0

        k = 0
        while k <= 10:
            k += 1
            # Input data is y/x -> needs switch
            result = lowess(numpy.array(data2), numpy.array(data1), delta=delta, frac=frac, it=10)

            if any( [math.isnan(r[1]) for r in result] ):
                print ("WARNING: lowess returned NA data points! We are trying to fix it")
                delta = delta * k
                result = lowess(numpy.array(data2), numpy.array(data1), delta=delta, frac=frac, it=10)
                frac = 1.0
            else:
                break

        return [ r[0] for r in result], [r[1] for r in result]

class LowessSmoothingCyLowess(LowessSmoothingBase):
    """Smoothing using Lowess smoother and then interpolate on the result
    """

    def __init__(self):
        pass

    def _initialize(self, data1, data2):
        try:
            import cylowess
            lowess = cylowess.lowess
        except ImportError:
            print("===================================")
            print("Cannot import the module lowess from 'cylowess', \nplease install the cylowess package according to http://slendermeans.org/lowess-speed.html (see also README)")
            print("===================================")

        delta = (max(data1) - min(data1)) * 0.01
        # Input data is y/x -> needs switch
        result = lowess(numpy.array(data2), numpy.array(data1), delta=delta, frac=0.1, it=10)
        return [ r[0] for r in result], [r[1] for r in result]

class UnivarSplineNoCV:
    """Smoothing of 2D data using a Python spline (no crossvalidation).

    Will use UnivariateSpline internally, it seems to have a tendency to
    overfit.
    """

    def __init__(self):
        pass

    def initialize(self, data1, data2):
        from scipy.interpolate import UnivariateSpline

        data1s, data2s = zip(*sorted(zip(data1, data2)))

        self.sp = UnivariateSpline(data1s, data2s)

    def predict(self, xhat):
        if len(xhat) == 1:
            return [ self.sp(xhat) ]
        else:
            return list(self.sp(xhat))

class UnivarSplineCV:
    """Smoothing of 2D data using a Python spline (using crossvalidation to determine smoothing parameters).

    Will use UnivariateSpline internally, setting the scipy smoothing parameter
    optimally "s" using crossvalidation  with part of the data (usually 25/75
    split). This prevents overfit to the data.
    """

    def __init__(self):
        pass

    def initialize(self, data1, data2, frac_training_data = 0.75, max_iter = 100, s_iter_decrease = 0.75, verb=False):
        from scipy.interpolate import UnivariateSpline

        if verb: 
            print(" --------------------" )
        
        # Random subsetting of parts of the data
        train_idx = random.sample(range(len(data1)), int(len(data1)*frac_training_data) )
        i = 0
        train_data1 = []
        train_data2 = []
        test_data1 = []
        test_data2 = []
        for d1,d2 in zip(data1, data2):
            if i in train_idx:
                train_data1.append(data1[i])
                train_data2.append(data2[i])
            else:
                test_data1.append(data1[i])
                test_data2.append(data2[i])
            i += 1

        # Sorted data points
        data1s, data2s = zip(*sorted(zip(data1, data2)))
        test_data1s, test_data2s = zip(*sorted(zip(test_data1, test_data2)))
        train_data1s, train_data2s = zip(*sorted(zip(train_data1, train_data2)))

        # Use initial linear Smoothing to find good smoothing parameter s
        smlin = SmoothingLinear()
        smlin.initialize(data2, data1)
        data2_lin_aligned = smlin.predict(data2)
        stdev_lin = numpy.std(numpy.array(data1) - numpy.array(data2_lin_aligned))
        linear_error = stdev_lin*stdev_lin

        # Perform initial spline approximation
        self.s = linear_error * len(train_data1s)
        self.sp = UnivariateSpline(train_data1s, train_data2s, k=3, s=self.s)

        # Apply spline approximation to the testdata
        test_data1_aligned = self.sp(test_data1)
        test_stdev = numpy.std(numpy.array(test_data2) - numpy.array(test_data1_aligned))
        if verb:
            test_median = numpy.median(numpy.array(test_data2) - numpy.array(test_data1_aligned))
            train_data1_aligned = self.sp(train_data1)
            tr_stdev = numpy.std(numpy.array(train_data2) - numpy.array(train_data1_aligned))
            tr_median = numpy.median(numpy.array(train_data2) - numpy.array(train_data1_aligned))
            print("  Lin:Computed stdev", stdev_lin)
            print("  Train Computed stdev", tr_stdev, "and median", tr_median)
            print("  Test Computed stdev", test_stdev, "and median", test_median)

        stdev_prev = test_stdev
        s_prev = self.s
        s_iter = self.s
        myIter = 0
        for i in range(max_iter):
            s_iter = s_iter * s_iter_decrease
            self.sp = UnivariateSpline(train_data1s, train_data2s, k=3, s=s_iter)
            test_data1_aligned = self.sp(test_data1)
            stdev = numpy.std(numpy.array(test_data2) - numpy.array(test_data1_aligned))
            if verb:
                print(" == Iter", s_iter, "\tstdev",  numpy.std(numpy.array(test_data2) - numpy.array(test_data1_aligned)))

            # Stop if stdev does not improve significantly any more
            #if stdev_prev - stdev < 0 or (i > 5 and (stdev_prev - stdev < 0.5)):
            if stdev_prev - stdev < 0:
                break

            stdev_prev = stdev
            s_prev = s_iter
            
        if verb:
            print(" == Done ", s_prev)

        # Final spline
        self.s = s_prev
        self.sp = UnivariateSpline(data1s, data2s, k=3, s=self.s)

    def predict(self, xhat):
        if len(xhat) == 1:
            return [ self.sp(xhat) ]
        else:
            return list(self.sp(xhat))

class SmoothingEarth:
    """Class for MARS type smoothing based on pyearth

    Get it at https://github.com/jcrudy/py-earth/
    """

    def __init__(self):
        pass

    def initialize(self, data1, data2):
        try:
            import pyearth
        except ImportError:
            print("===================================")
            print("Cannot import pyearth, \nplease install it from https://github.com/jcrudy/py-earth/")
            print("===================================")

        self.model = pyearth.Earth()
        X = numpy.array(data1)
        Y = numpy.array(data2)
        self.model.fit(X, Y)

    def predict(self, xhat):
        return list(self.model.predict(xhat))

class SmoothingLinear:
    """Class for linear transformation
    """

    def __init__(self):
        pass

    def initialize(self, data1, data2):
        # data1 is the predictor (e.g. the input) -> x
        # data2 is the response (e.g. what we want to predict) -> y
        A = numpy.array([ numpy.array(data1), numpy.ones(len(data1))])
        self.w = numpy.linalg.lstsq(A.T,numpy.array(data2))[0] # obtaining the parameters

    def predict(self, xhat):
        xhat_np = numpy.array(xhat)
        predicted_result = self.w[0]*xhat_np+self.w[1] # regression line
        return list(predicted_result)

class SmoothingInterpolation:
    """Class for interpolation transformation
    """

    def __init__(self):
        pass

    def initialize(self, data1, data2):
        # data1 is the predictor (e.g. the input) -> x
        # data2 is the response (e.g. what we want to predict) -> y
        self.use_cy = True

        data1s, data2s = zip(*sorted(zip(data1, data2)))

        try:
            from msproteomicstoolslib.cython._optimized import CyLinearInterpolateWrapper
            self.f = CyLinearInterpolateWrapper(data1s, data2s, 0.0)
        except ImportError:
            print("WARNING: cannot import CyLinearInterpolateWrapper, will use Python version (slower).")
            from scipy.interpolate import interp1d
            self.f = interp1d(data1s, data2s)
            self.use_cy = False

        # Also prepare linear transformation
        self.linear_sm = SmoothingLinear()
        self.linear_sm.initialize(data1, data2)

    def getLWP(self):
        if self.use_cy:
            return self.f
        else:
            raise Exception("Cannot return CyLinearInterpolateWrapper wrapper!")

    def predict(self, xhat):
        try:
            if self.use_cy:
                predicted_result = self.f.predict(list(xhat)) # interpolation fxn
                return predicted_result
            else:
                predicted_result = self.f(xhat) # interpolation fxn
        except ValueError:
            # outside bound, use linear
            return self.linear_sm.predict(xhat)

        # If the input is a list, we need to check every value in the list
        # -> fix those that are "NA"
        if any( [math.isnan(pp) for pp in predicted_result] ):
            
            for i, (qq, pp) in enumerate(zip(xhat, predicted_result)):
                if (math.isnan(pp)):
                    predicted_result[ i ] = self.linear_sm.predict([ qq ])[0]


        return list(predicted_result)

class LocalKernel:
    """Base class for local kernel smoothing
    """

    def initialize(self, data1, data2):
        # data1 is the predictor (e.g. the input) -> x
        # data2 is the response (e.g. what we want to predict) -> y
        data1, data2 = zip(*sorted(zip(data1, data2)))
        self.data1 = numpy.array(data1)
        self.data2 = numpy.array(data2)

        if self.removeOutliers:
            pass

    def _getLocalDatapoints(self, data1, data2, topN_, max_diff, xhat):
            """ Return all datapoints that are within max_diff of xhat

            If there are less than 2 * topN datapoints within max_diff of xhat,
            returns the 2 * topN datpoints around xhat.
            """

            topN = int(topN_) # ensure that topN is integer

            if len(data1) < 2*topN:
                return data1, data2

            # This lower bound will actually get the element that is just larger
            # than the search parameter
            lb = abs(numpy.searchsorted(data1, xhat)) - 1
            if lb - topN < 0:
                lb = topN
            if lb >= len(data1):
                lb = len(data1)-topN

            source_d = []
            target_d = []

            # Walk to the left (down)
            it = lb
            while it >= 0:
                if abs(data1[it] - xhat) < max_diff:
                    source_d.append(data1[it])
                    target_d.append(data2[it])
                else: 
                    break
                it -= 1

            # Walk to the right (up)
            it = lb + 1
            while it < len(data1):
                if abs(data1[it] - xhat) < max_diff:
                    source_d.append(data1[it])
                    target_d.append(data2[it])
                else: 
                    break
                it += 1

            # Check if we have enough datapoints
            if len(source_d) < 2*topN:
                lower = lb-topN if lb-topN >= 0 else 0
                upper = lb+topN if lb+topN < len(data1) else len(data1)
                return numpy.asarray(data1[lower:upper]), numpy.asarray(data2[lower:upper])
            else:
                return numpy.asarray(source_d), numpy.asarray(target_d)

class WeightedNearestNeighbour(LocalKernel):
    """Class for weighted interpolation using local linear differences

    This function uses the weighted mean of the k nearest neighbors to
    calculate the transformation.  This method may be affected by single
    outlier close to the transformation point.

    Each neighboring point is given a weight equal to  ::

		   1
	--------------------------
	  abs( distance ) ** exp 


    up to a minimal distance min_diff after which the weight cannot increase any more.
    """

    def __init__(self, topN, max_diff, min_diff, removeOutliers, exponent=1.0):
        """ Initialize WNN

        Args:
            topN(integer): how many datapoints should be included for nearest neighbor
            max_diff(float): maximal difference in x to be included for nearest neighbor
            min_diff(float): minimal difference in x to still give proportional weight
            removeOutliers(bool): no effect
            exponent(float): exponent to be used
        """
        assert topN is not None 

        self.topN = topN
        self.max_diff = max_diff
        self.min_diff = min_diff
        self.removeOutliers = removeOutliers
        self.EXP = exponent

    def predict(self, xhat):

        res = []
        for xhat_ in xhat:

            source_d, target_d = self._getLocalDatapoints(
                self.data1, self.data2, self.topN, self.max_diff, xhat_)

            # Transform target data:
            #   Compute a difference array from the source and apply it to the target
            #   (local linear differences)
            source_d_diff = source_d - xhat_
            target_data_transf = target_d - source_d_diff

            # Use transformed target data to compute expected RT in target domain (weighted average)
            # EXP = 0.65 # produces slightly better results
            EXP = self.EXP
            weights = numpy.clip( numpy.abs(source_d_diff), self.min_diff, numpy.inf)
            weights = 1. / weights**EXP
            expected_targ = numpy.average(target_data_transf, weights=weights)

            # Compute a measurement of dispersion, standard deviation
            self.last_dispersion = numpy.std(target_data_transf)

            res.append( expected_targ )

        return res

class SmoothLLDMedian(LocalKernel):
    """Class for local median interpolation using local linear differences

    This function uses the median of the k nearest neighbors to calculate the
    transformation.  This is robust, unweighted method as a single outlier will
    not substantially affect the result.

    This method assumes that the data is locally smooth and linear
    """

    def __init__(self, topN, max_diff, min_diff, removeOutliers):
        assert topN is not None 

        self.topN = topN
        self.max_diff = max_diff
        self.min_diff = min_diff
        self.removeOutliers = removeOutliers

    def predict(self, xhat):

        def mad(data, axis=None):
            """Median absolute deviation (MAD) is a robust estimator of variation.

            http://en.wikipedia.org/wiki/Robust_measures_of_scale#IQR_and_MAD
            """
            return median(absolute(data - median(data, axis)), axis)

        res = []
        for xhat_ in xhat:

            source_d, target_d = self._getLocalDatapoints(self.data1, self.data2, self.topN, self.max_diff, xhat_)

            # Transform target data:
            #   Compute a difference array from the source and apply it to the target
            #   (local linear differences)
            source_d_diff = [s - xhat_ for s in source_d]
            target_data_transf = [t - s for t,s in zip(target_d, source_d_diff)]

            # Use transformed target data to compute expected RT in target domain (median)
            expected_targ = numpy.median(target_data_transf)

            # Compute a (robust) measurement of dispersion, median absolute deviation (MAD)
            self.last_dispersion = mad(target_data_transf)

            res.append( expected_targ )

        return res

