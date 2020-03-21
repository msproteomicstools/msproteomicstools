#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta ChromatogramSmoothing.py>
Copyright (C) 2020 Shubham Gupta

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

--------------------------------------------------------------------------
$Maintainer: Shubham Gupta$
$Authors: Shubham Gupta$
--------------------------------------------------------------------------
"""
import numpy as np
import math
from msproteomicstoolslib.math.Smoothing import LowessSmoothingBase
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

def getXIC_SmoothingObj(smoother, kernelLen = 13, polyOrd = 4):
    if smoother == "sgolay":
        return SgolaySmooth(kernelLen, polyOrd)
    elif smoother == "gaussian":
        return GaussianSmooth(kernelLen)
    elif smoother == "loess":
        return LoessSmooth(kernelLen, polyOrd)
    else:
        raise Exception("Unknown chromatogram smoothing method: " + smoother)

class LoessSmooth(LowessSmoothingBase):
    """Smoothing chromatogram using Lowess smoother and then interpolate on the result
    """

    def __init__(self, kernelLen, polyOrd):
        self.kernelLen = kernelLen
        self.polyOrd = polyOrd

    def _initialize(self, data1, data2):
        try:
            import statsmodels.api as sm
            lowess = sm.nonparametric.lowess
        except ImportError:
            print("===================================")
            print("Cannot import the module lowess from 'statsmodels', \nplease install the Python package 'statsmodels'")
            print("===================================")
        
        if len(data2) < 100:
            frac = 1.0
        else:
            frac = self.kernelLen / len(data1)

        k = 0
        while k <= 10:
            k += 1
            # Input data is y/x -> needs switch
            result = lowess(np.array(data2), np.array(data1), delta=0.0, frac=frac, it=0)

            if any( [math.isnan(r[1]) for r in result] ):
                print ("WARNING: lowess returned NA data points! We are trying to fix it")
                result = lowess(np.array(data2), np.array(data1), delta=0.0, frac=frac, it=10)
                frac = 1.0
            else:
                break

        return [ r[0] for r in result], [r[1] for r in result]
    
    def smooth(self, xhat, yhat):
        return xhat, np.array(self.predict(xhat))

class SgolaySmooth:
    """Savitzky-Golay smoothing. It preserves the peak shape"""
    def __init__(self, kernelLen, polyOrd):
        if (kernelLen % 2) == 0:
            raise Exception("kernel length must be odd for Savitzky- Golay smoothing.")
        if not (kernelLen > polyOrd):
            raise Exception("polyOrd must be less than the kernel length for Savitzky- Golay smoothing.")
        self.window_length = kernelLen
        self.polyorder = polyOrd

    def initialize(self, data1, data2):
        pass

    def smooth(self, xhat, yhat):
        return xhat, savgol_filter(np.array(yhat), self.window_length, self.polyorder)

class GaussianSmooth:
    """
    Gaussian smoothing. The kernelLen indicates Full Width Half Maximum (FWHM) of gaussian kernel.
    The standard deviation is multiplied with 0.3706505 to be similar to R code
    https://github.com/shubham1637/DIAlignR/R/chromatogram_smooth.R

    """
    def __init__(self, kernelLen):
        # Converting FWHM to standard deviation.
        sigma = kernelLen / np.sqrt(8 * np.log(2))
        self.sigma = sigma*0.3706505 # Scaled to have quartiles at +/- 0.25*sigma

    def initialize(self, data1, data2):
        pass

    def smooth(self, xhat, yhat):
        return xhat, gaussian_filter1d(np.array(yhat), sigma = self.sigma, mode = 'constant', cval = 0.0)
