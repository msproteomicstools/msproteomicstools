#!/usr/bin/env python
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
$Maintainer: Hannes Roest $
$Authors: Hannes Roest $
--------------------------------------------------------------------------
"""
from __future__ import print_function


import sys, re
import numpy
import pyopenms

if len(sys.argv) < 4:
  print("A small program that will filter chromatogram mzML file by their native id") 
  print("Usage: python filterChrom.py infile.mzML outfile.mzML chrom_id inverse") 
  print("Note: chrom_id can be a regex; using True as 'inverse' argument allows to do inverse matching")
  sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]
filter_criteria = sys.argv[3]
inverse = False
if len(sys.argv) > 4:
    inverse = bool(sys.argv[4])

print("Will filter with criteria ", filter_criteria, "(inverse: %s)" % inverse)

chroms_out = []
try:
    # PyMzML path
    import pymzml
    exp2 = pyopenms.MSExperiment()
    run = pymzml.run.Reader(infile, build_index_from_scratch=True)
    for key in list(run.info['offsets'].keys()):
        if key == "indexList": continue
        if key == "TIC": continue
        if (inverse and not re.search(filter_criteria, key)) \
           or (not inverse and re.search(filter_criteria, key)):
            c = pyopenms.MSChromatogram()
            c.setNativeID(str(key))
            pr = pyopenms.Precursor()
            chrom = run[key]
            #try: 
            pr.setMZ( chrom["precursors"][0]["mz"] )
            #except Exception: pass
            c.setPrecursor(pr)
            timea = numpy.array( chrom.time , dtype=numpy.float32)
            inta = numpy.array( chrom.i , dtype=numpy.float32)
            peaks = numpy.ndarray(shape=(len(timea), 2), dtype=numpy.float32)
            peaks[:,0] = timea
            peaks[:,1] = inta
            c.set_peaks(peaks)
            chroms_out.append(c)

except ImportError:

    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(infile, exp)
    exp2 = exp
    exp2.clear(False)
    chroms = exp2.getChromatograms()
    for c in chroms:
        if (inverse and not re.search(filter_criteria, key)) \
           or (not inverse and re.search(filter_criteria, key)):
            chroms_out.append(c)


# Sort chromatograms and store again
print("Retrieved", len(chroms_out), "chromatograms.")
chroms_out.sort(key=lambda x: x.getNativeID())
exp2.setChromatograms(chroms_out)
pyopenms.MzMLFile().store(outfile, exp2)
