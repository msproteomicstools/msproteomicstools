#!/usr/bin/python
# -*- coding: utf-8 -*-
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


class SqlDataAccess(object):

    def __init__(self, filename):
        import sqlite3
        self.conn = sqlite3.connect(filename)
        self.c = self.conn.cursor()

    def getDataForChromatograms(self, ids):
        """
        Get data from multiple chromatograms chromatogram

        - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        - data_type is one of 0 = mz, 1 = int, 2 = rt
        - data contains the raw (blob) data for a single data array
        """

        stmt ="SELECT CHROMATOGRAM_ID, COMPRESSION, DATA_TYPE, DATA FROM DATA WHERE CHROMATOGRAM_ID IN ("
        for myid in ids:
            stmt += str(myid) + ","
        stmt = stmt[:-1]
        stmt += ")"

        data = [row for row in self.c.execute(stmt)]
        tmpres = self._returnDataForChromatogram(data)

        res = []
        for myid in ids:
            res.append( tmpres[myid] )
        return res

    def getDataForChromatogram(self, myid):
        """
        Get data from a single chromatogram

        - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        - data_type is one of 0 = mz, 1 = int, 2 = rt
        - data contains the raw (blob) data for a single data array
        """

        data = [row for row in self.c.execute("SELECT CHROMATOGRAM_ID, COMPRESSION, DATA_TYPE, DATA FROM DATA WHERE CHROMATOGRAM_ID = %s" % myid )]
        return self._returnDataForChromatogram(data).values()[0]

    def getDataForChromatogramFromNativeId(self, native_id):
        """
        Get data from a single chromatogram

        - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        - data_type is one of 0 = mz, 1 = int, 2 = rt
        - data contains the raw (blob) data for a single data array
        """

        data = [row for row in self.c.execute("SELECT CHROMATOGRAM_ID, COMPRESSION, DATA_TYPE, DATA FROM DATA INNER JOIN CHROMATOGRAM ON CHROMATOGRAM.ID = CHROMATOGRAM_ID WHERE NATIVE_ID = %s" % native_id )]
        return self._returnDataForChromatogram(data).values()[0]

    def _returnDataForChromatogram(self, data):
        import PyMSNumpress
        import zlib

        # prepare result
        chr_ids = set([chr_id for chr_id, compr, data_type, d in data] )
        res = { chr_id : [None, None] for chr_id in chr_ids }

        rt_array = []
        intensity_array = []
        for chr_id, compr, data_type, d in data:
            result = []
            if compr == 5:
                tmp = [ord(q) for q in zlib.decompress(d)]
                PyMSNumpress.decodeLinear(tmp, result)
            if compr == 6:
                tmp = [ord(q) for q in zlib.decompress(d)]
                PyMSNumpress.decodeSlof(tmp, result)

            if data_type == 1:
                res[chr_id][1] = result
            elif data_type == 2:
                res[chr_id][0] = result
            else:
                raise Exception("Only expected RT or Intensity data for chromatogram")

        return res
