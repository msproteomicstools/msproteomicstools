#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta smooth_chromatograms.py>
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

import msproteomicstoolslib.math.ChromatogramSmoothing as xic_smoothing

class chromSmoother():
    """
    Use the ChromatogramSmoothing part of msproteomicstoolslib to smooth a chromatogram-group.

    >>> chrom_smoother = chromSmoother()
    >>> XIC_Group_sm = chrom_smoother.rt_smoothXICs(XIC_Group)
    """
    def __init__(self, smoother="sgolay", kernelLen=11, polyOrd = 4):
        """
        Initializes ChromatogramSmoothing with either sgolay, gaussian or loess.
        kernelLen defines the numbers of neighboring-points to be used by the smoothing function.
        polyOrd defines the order of the polynomial fit used by the smoothing function.
        """
        self.sm = xic_smoothing.getXIC_SmoothingObj(smoother, kernelLen, polyOrd)

    def smoothXICs(self, XIC_Group):
        """
        Smooth each element of XIC_Group.

        >>> import math
        >>> XIC_Group = [(np.array([i for i in range(15)]), np.array([1+math.sin(i) for i in range(15)]))]
        >>> XIG_Group_sm = chrom_smoother.smoothXICs(XIC_Group)
        """
        XIC_Group_sm = []
        for XIC in XIC_Group:
            self.sm.initialize(XIC[0], XIC[1])
            XIC_sm = self.sm.smooth(XIC[0], XIC[1])
            XIC_Group_sm.append(XIC_sm)
        return XIC_Group_sm
