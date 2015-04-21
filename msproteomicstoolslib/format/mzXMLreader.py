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
from xml.etree.cElementTree import iterparse
import xml.etree.cElementTree as etree
import numpy #for the "floor" fxn
import sys
try:
    from StringIO import StringIO
except ImportError:
    # Python 3
    from io import StringIO


from ..data_structures import DDB


"""
Doc :
    A small class to read a subset of the mzXML file.
"""


class mzXMLReader:

    def __init__(self, filename, init = False):
        self.file = filename
        self.scan_ms_bins = None
        self.scan_rt_bins = None
        self.index        = None
        if init:
            self.parse_scans(ms2Only = True, createHashes = True)

    def parse_header(self):
        try:
            #about 1000x faster, it can distinguish non mzXML files
            self._parse_header_fast()
        except SyntaxError as e:
            #print e.args
            raise Exception('Parsing failed, probably not an mzXML file')

    def _parse_header_fast(self):
        f = open( self.file, 'r')
        lines = []
        while True:
            line = f.readline()
            if line == None or line.find( "<scan"  ) != -1: break
            lines.append( line )

        string = "".join( lines )
        string += "\n</msRun>"
        string += "\n</mzXML>"
        source = StringIO( string )
        context = iterparse(source, events=("start", "end", "start-ns"))
        self._parse_header_context(context)

    def _parse_header_context(self, context):
        doneI = False
        doneP = False
        for event, elem in context:
            try:
                if elem.tag[ -12: ] == 'msInstrument':
                    doneI = True
                    for child in elem.getchildren():
                        if child.get('category') == 'msManufacturer':   
                            self.manu = child.get('value')
                        if child.get('category') == 'msModel':          
                            self.model = child.get('value')
                        if child.get('category') == 'msIonisation':     
                            self.ioni = child.get('value')
                        if child.get('category') == 'msMassAnalyzer':   
                            self.analy = child.get('value')
                        if child.get('category') == 'msDetector':       
                            self.detect = child.get('value')
                        if child.tag[-8:] == 'software':    
                            self.instrumentSoftwareType = child.get('type')
                            self.instrumentSoftwareName = child.get('name')
                            self.instrumentSoftwareVersion = child.get(
                                    'version')
                if elem.tag[ -14: ] == 'dataProcessing':
                    doneP = True
                    for child in elem.getchildren():
                        if child.tag[-8:] == 'software':    
                            self.analyseSoftwareType = child.get('type')
                            self.analyseSoftwareName = child.get('name')
                            self.analyseSoftwareVersion = child.get( 
                                    'version')
            except Exception: pass
            if doneP and doneI: break

    def _parse_header(self):
        #f = open( self.file, 'r')
        #source = open( rootdir + testdir + testfile, 'r')
        source = open( self.file, 'r')
        context = iterparse(source, events=("start", "end", "start-ns"))
        self._parse_header_context()

    def parse_scans(self, ms2Only = True, createHashes = False, readPeaks=False):
        self.index = {}
        self.scan_number_hash = {}
        scans = []
        fail = []
        source = open( self.file, 'r')
        context = iterparse(source, events=("start", "end"))
        for event, elem in context:
            if event == "end" and elem.tag[ -4: ] == 'scan': 
                s = Scan( elem)
                # print  "scan ", s, s.msLevel, ms2Only, s.msLevel == 2 and ms2Only
                if s.msLevel == 2 and len(elem) ==0: fail.append( elem )
                elif s.msLevel == 2 and ms2Only: 
                    if readPeaks: s.read_peaks( elem )
                    scans.append( s )
                    self.scan_number_hash[ s.scan_number ] = s
                else:
                    # if readPeaks: s.read_peaks( elem )
                    scans.append( s )
                    self.scan_number_hash[ s.scan_number ] = s
            if event == "end" and elem.tag.endswith("index"):
                for child in elem.getiterator():
                    if child.tag.endswith( "offset" ):
                        self.index[ int(child.get( "id" )) ] = int(child.text)

        if createHashes: self.create_ms_rt_hashes( scans )
        return scans

    def create_ms_rt_hashes(self, scans):
        self.scan_ms_bins = {}
        self.scan_rt_bins = {}
        for s in scans:
            prec_floor = numpy.floor( s.precursorMZ )
            rt_floor = numpy.floor( s.retTime )
            if prec_floor in self.scan_ms_bins: 
                self.scan_ms_bins[ prec_floor ].append( s)
            else: self.scan_ms_bins[ prec_floor ] = [ s ]
            if rt_floor in self.scan_rt_bins:
                self.scan_rt_bins[ rt_floor ].append( s )
            else: self.scan_rt_bins[ rt_floor ] = [s]

    def read_scan(self, id, readPeaks = False):
        if self.index is None: self.parse_scans(ms2Only = True, 
                createHashes = True)
        position = self.index[id]
        source = open( self.file, 'r')
        source.seek( position ) 
        #read in, line based
        lines = []
        for line in source.readlines():
            lines.append( line )
            if line.strip() == "</scan>": break
        #cast to string and feed to parses StringIO based :-)
        string = "".join( lines )
        source = StringIO( string )
        context = iterparse(source, events=("start", "end") )
        for event, elem in context:
            if event == "end" and elem.tag[ -4: ] == 'scan': 
                myscan = Scan( elem)
                if readPeaks: myscan.read_peaks( elem )
        return myscan

    def find_corresponding_scan( self, searchHit,
            readScan = False, readPeaks = False):

        #raise ValueError("Implement search based on scan number")
        scan = self.scan_number_hash[ searchHit.scan_number ]
        scan.searchHit = searchHit
        searchHit.scan = scan
        return  scan

        if self.scan_rt_bins is None: 
            raise ValueError('you need to run create_ms_rt_hashes() first')

        rt = searchHit.get_retention_time()
        m  = searchHit.get_precursor_mass() 
        #find all scans that are in the same RT window
        retWindow = self.scan_rt_bins[ numpy.floor(rt) ] 
        for scan in retWindow:
            if abs(scan.retTime - rt) < 0.05: break

        #assert we have the right one
        assert abs(scan.precursorMZ - m ) < 0.1

        if readScan: scan = self.read_scan( scan.id, readPeaks)
        return scan

class Scan:
    def __init__(self, elem, is_mrm = False):
        self.id        = int(   elem.get("num"))
        self.msLevel   = int(   elem.get("msLevel"))
        self.peaksC    = int(   elem.get("peaksCount"))
        self.retTime   = round( float( elem.get("retentionTime")[2:-1] ), 2)
        self.basePeak  = float( elem.get("basePeakMz"))
        self.basePeakI = float( elem.get("basePeakIntensity"))
        self._peaks    = None
        self.precursorMZ = -1
        self.mrm       = is_mrm
        for child in elem.getchildren():
            if child.tag[-11 :] == 'precursorMz':
                #try:
                self.precursorIntensity = float(child.get('precursorIntensity'))
                self.activationMethod   = child.get('activationMethod')
                self.precursorMZ        = float( child.text )
                try:
                    self.precursorCharge    = float(child.get('precursorCharge'))
                except Exception: 
                    self.mrm = True

    def read_peaks(self, elem):
        self._peaks = []
        string = ''
        for child in elem.getchildren():
            if child.tag.endswith( 'peaks' ):
                self.precision = child.get( 'precision' )
                self.byteOrder = child.get( 'byteOrder' )
                self.pairOrder = child.get( 'pairOrder' )
                if not self.pairOrder:
                    self.pairOrder = child.get('contentType')
                string = child.text
        assert (self.byteOrder == 'network' and self.precision == '32' and
        self.pairOrder == 'm/z-int')
        decoder = mzXML64coder()
        peaks_decoded = decoder.decode(string)
        self._peaks = [ Peak( p[0] , p[1] ) for p in peaks_decoded ]

    @property
    def peaks(self):
        if self._peaks is None: raise ValueError("""You have to read in\
 the peaks first before accessing. Use the read_scan function on the reader
 object.""")
        else: return self._peaks

    @property 
    def scan_number(self):
        return self.id

    def normalize_peaks(self, factor = 10000):
        max = self.max_peak()
        for p in self.peaks:
            p.int = p.int * 10000/ max.int

    def max_peak(self):
        max = self.peaks[0]
        for p in self.peaks:
            if p.int > max.int: max = p
        return max

    def annotate_peaks(self, myrange=range(1,4), do_phospho = False, do_ax=False):
        R = silver.Residues.Residues('mono')
        peptide = DDB.Peptide()
        peptide.set_sequence( self.searchHit.modified_peptide )
        pepLen = len( peptide )
        mz_tol = 0.7
        k = 0
        for p in self.peaks:
            tests = 0
            p.min = 99
            S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
            S.ion_charge = self.searchHit.get_assumed_charge()
            S.construct_from_peptide( peptide.get_modified_sequence('SEQUEST'), 
                                     R.residues, R.res_pairs)
            for ch in myrange:
                S.ion_charge = ch
                do_phospho = peptide.has_phospho()
                annotate_peak( p, ch, S, R, tests, pepLen, do_phospho)
            if abs(p.min) > mz_tol: p.deannotate()
            else: k+= 1
            p.annotation(), p.min, p.mz, p.int

    def print_peaks_above(self, cut_at):
        #print annotation
        cut_at = 10
        self.peaks.sort( lambda x,y: cmp( x.int, y.int) )
        #scan.peaks.sort( lambda x,y: cmp( x.mz, y.mz) )
        maxInt = max(  [ p.int for p in self.peaks] )
        for p in self.peaks:
            if p.int > maxInt / cut_at: ( [str(p.annotation('spectrast')),
                str(p.mz), str(p.int)  ] )

def annotate_peak(p, ch, S, R, tests, pepLen, do_phospho=False, do_ax=False):
        #neutral loss from the parent ion
        if do_phospho:
            pred = (S.peptide_mass - S.mass_H1PO3) / ch
            if abs( p.mz - pred ) < abs(p.min):
                p.min = p.mz - pred
                p.annotate( 'pn', ch, 0, S.mass_H1PO3,  'HPO3' )
            pred = (S.peptide_mass - S.mass_H3PO4) / ch
            if abs( p.mz - pred ) < abs(p.min):
                p.min = p.mz - pred
                p.annotate( 'pn', ch, 0, S.mass_H3PO4, 'H3PO4' )
        #pred = (S.peptide_mass - R.mass_OH - R.mass_H) / ch
        #print pred
        #if abs( p.mz - pred ) < abs(p.min):
        #    p.min = p.mz - pred
        #    p.annotate( 'pn', ch, 0, S.mass_OH + S.mass_H , 'H2O' )

        #y, b, a, x series
        charged_y_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.y_series ] 
        charged_b_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.b_series ] 
        charged_a_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.a_series ] 
        charged_x_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.x_series ] 
        charged_n_series =  [ ( pred + (ch -1)*S.mass_H)/ch for pred in S.neutral_losses ] 

        for i, pred in enumerate(charged_y_series):
            tests += len( charged_y_series )
            if abs( p.mz - pred ) < abs(p.min):
                p.min = p.mz - pred
                p.annotate( 'y', ch, pepLen - 1 - i)
        for i, pred in enumerate(charged_b_series):
            tests += len( charged_y_series )
            #adj_pred = (pred + (ch -1)*S.mass_H ) / ch
            if abs( p.mz - pred ) < abs(p.min):
                p.min = p.mz - pred
                p.annotate( 'b', ch, i + 1)

        #lets also look at y and b neutral losses of phospho
        if do_phospho:
            tests += 2 + 2 * len( S.y_series )
            for i, pred in enumerate(S.y_series):
                pred = (pred + (ch -1)*S.mass_H - S.mass_H1PO3) / ch
                if abs( p.mz - pred ) < abs(p.min):
                    p.min = p.mz - pred
                    p.annotate( 'yn', ch, pepLen - 1 -i, S.mass_H1PO3, 'HPO3')
            for i, pred in enumerate(S.y_series):
                pred = (pred + (ch -1)*S.mass_H - S.mass_H3PO4) / ch
                if abs( p.mz - pred ) < abs(p.min):
                    p.min = p.mz - pred
                    p.annotate( 'yn', ch, pepLen - 1 -i, S.mass_H3PO4, 'H3PO4')
            tests += 2 + 2 * len( S.b_series )
            for i, pred in enumerate(S.b_series):
                pred = (pred + (ch -1)*S.mass_H - S.mass_H1PO3) / ch
                if abs( p.mz - pred ) < abs(p.min):
                    p.min = p.mz - pred
                    p.annotate( 'bn', ch, i +1, S.mass_H1PO3, 'HPO3')
            for i, pred in enumerate(S.b_series):
                pred
                pred = (pred + (ch -1)*S.mass_H - S.mass_H3PO4) / ch
                i, pred
                if abs( p.mz - pred ) < abs(p.min):
                    p.min = p.mz - pred
                    p.annotate( 'bn', ch, i + 1, S.mass_H3PO4, 'H3PO4')

        #the a and x ion series
        if do_ax:
            for i, pred in enumerate(charged_a_series):
                tests += len( charged_y_series )
                if abs( p.mz - pred ) < abs(p.min):
                    p.min = p.mz - pred
                    p.annotate( 'a', ch, i + 1)
            for i, pred in enumerate(charged_x_series):
                tests += len( charged_y_series )
                if abs( p.mz - pred ) < abs(p.min):
                    p.min = p.mz - pred
                    p.annotate( 'x', ch, pepLen - i - 1 )

        #neutral loss of H20 or NH3
        for i, pred in enumerate(charged_n_series):
            tests += 2* len( charged_y_series )
            if abs( p.mz - pred ) < abs(p.min):
                if i % 2 == 0:
                    p.min = p.mz - pred
                    p.annotate( 'bn', ch, i/2 + 1, S.mass_H2O, 'H2O')
                else:
                    p.min = p.mz - pred
                    p.annotate( 'yn', ch, pepLen - i/2 -1 , S.mass_NH3, 'NH3')

class Peak:
    def __init__(self, mz, int, peptide = None):
        self.mz = mz
        self.int = int
        self.annotated = False
        self.isotope = 0

    def annotate(self, type, charge, number, loss = 0, comment = ''):
        self.annotated = True
        self.type = type
        self.charge = charge
        self.number = number
        self.loss = loss
        self.comment = comment

    def annotate_isotope(self, peak, plus):
        self.annotated = True
        self.type =     peak.type
        self.charge =   peak.charge
        self.number =   peak.number
        self.loss =     peak.loss
        self.comment =  peak.comment
        self.isotope =  peak.isotope + plus

    def deannotate(self):
        self.annotated = False

    #@property
    def annotation(self, format = 'viewer'):
        if not self.annotated: return 'None'
        format = format.strip()
        ch = ''
        iso = ''
        if self.isotope != 0: iso = 'i'
        if format == 'spectrast':
            if self.charge > 1: ch = "^%s" % self.charge
            if self.type in ('bn' , 'yn'):
                return "%s%s-%s%s%s/%s" % (self.type[0], self.number, self.loss,
                        ch, iso, self.min)
            if self.type in ('pn'):
                return "%s-%s%s%s/%s" % (self.type[0], self.loss,
                        ch, iso, self.min)
            return "%s%s%s%s/%s" % (self.type, self.number, ch, iso, self.min)
        elif format == 'viewer':
            for i in range(self.charge): ch += '+'
            if self.type == 'bn' or self.type == 'yn':
                return "(%s%s%s - %s)" % (self.type[0], self.number, ch,
                        self.comment )
            if self.type == 'pn':
                return "([M+%sxH]%s - %s)" % (self.charge, ch, self.comment )
            return "%s%s%s" % (self.type, self.number, ch)
        else: raise ValueError("format %s not known" % format)

    def is_y(self):
        return self.annotated and self.type in ('y', 'yn')

    def is_b(self):
        return self.annotated and self.type in ('b', 'bn')

    def is_parent(self):
        return self.annotated and self.type in ('p', 'pn')


class mzXML64coder:
    def __init__(self):
        pass

    def decode(self, mystring ):
        import base64
        import struct

        # Python 3 handles strings differently 
        if (sys.version_info > (3, 0)):
            mynr = base64.standard_b64decode(bytes(mystring, 'utf-8'))
        else:
            mynr = base64.standard_b64decode(mystring)

        peaks_ints = []
        for i in range(0, len(mynr), 8):
            peak = struct.unpack('>f', mynr[i:i+4]) 
            inte = struct.unpack('>f', mynr[i+4:i+8]) 
            peaks_ints.append( tuple( [peak[0], inte[0] ]) )
        return peaks_ints

    #give a list of tuples: ( peak, intensity )
    def encode(self, peaks_ints ):
        import base64
        import struct
        str = ''
        for t in peaks_ints:
            peak = t[0]
            ints = t[1]
            pStr = struct.pack( '>f', peak )
            iStr = struct.pack( '>f', ints )
            str += pStr + iStr
        str = base64.standard_b64encode( str )
        return str

def TEST_mzXMLdecoder():
    input ='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'
    coder = mzXML64coder()
    res = coder.decode( input )
    res2 = coder.encode( res )
    assert res2 == input
    print("Tests for mzXML decoder passed.")

if __name__ == '__main__':
    TEST_mzXMLdecoder()

