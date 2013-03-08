"""
Copyright (c) 2010-2012, Kapteyn Astronomical Institute, University
of Groningen. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Kapteyn Astronomical Institute nor the names
      of its contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Chauvenet from http://www.astro.rug.nl/software/kapteyn/plot_directive/EXAMPLES/kmpfit_hubblefit.py
from scipy.special import erfc
import numpy

def chauvenet(x, y, mean=None, stdv=None):
   #-----------------------------------------------------------
   # Input:  NumPy arrays x, y that represent measured data
   #         A single value of a mean can be entered or a 
   #         sequence of means with the same length as 
   #         the arrays x and y. In the latter case, the 
   #         mean could be a model with best-fit parameters.
   # Output: It returns a boolean array as filter.
   #         The False values correspond to the array elements
   #         that should be excluded
   # 
   # First standardize the distances to the mean value
   # d = abs(y-mean)/stdv so that this distance is in terms
   # of the standard deviation.  
   # Then the  CDF of the normal distr. is given by 
   # phi = 1/2+1/2*erf(d/sqrt(2))
   # Note that we want the CDF from -inf to -d and from d to +inf.
   # Note also erf(-d) = -erf(d).
   # Then the threshold probability = 1-erf(d/sqrt(2))
   # Note, the complementary error function erfc(d) = 1-erf(d)
   # So the threshold probability pt = erfc(d/sqrt(2))
   # If d becomes bigger, this probability becomes smaller.
   # If this probability (to obtain a deviation from the mean)
   # becomes smaller than 1/(2N) than we reject the data point
   # as valid. In this function we return an array with booleans 
   # to set the accepted values.
   # 
   # use of filter:
   # xf = x[filter]; yf = y[filter]
   # xr = x[~filter]; yr = y[~filter]
   # xf, yf are cleaned versions of x and y and with the valid entries
   # xr, yr are the rejected values from array x and y
   #-----------------------------------------------------------
   if mean is None:
      mean = y.mean()           # Mean of incoming array y
   if stdv is None:
      stdv = y.std()            # Its standard deviation
   N = len(y)                   # Lenght of incoming arrays
   criterion = 1.0/(2*N)        # Chauvenet's criterion
   d = abs(y-mean)/stdv         # Distance of a value to mean in stdv's
   d /= 2.0**0.5                # The left and right tail threshold values
   prob = erfc(d)               # Area normal dist.    
   filter = prob >= criterion   # The 'accept' filter array with booleans
   return filter                # Use boolean array outside this function

