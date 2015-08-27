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

from xml.etree.cElementTree import iterparse
import xml.etree.cElementTree as etree


"""
Doc :
    A small class to read a subset of the pep.xml file.
"""

class pepXMLReader:

    def __init__(self, filename):
        self.file = filename

    def parse_for_FDR(self, threshold = 0):
        try:
            return self._parse_for_FDR(threshold)
        except SyntaxError as e:
            raise Exception('Parsing failed')

    def parse_all(self):
        return self._parse_for_FDR(-1, True)

    def _parse_for_FDR(self, threshold, getAll=False):
        source = open( self.file, 'r')
        context = iterparse(source)
        search_hits = []
        for event, elem in context:
            if elem.tag.find("spectrum_query") != -1:
                specQuery = Spectrum_Query( elem )
                if len( elem.getchildren() ) >0:
                    temp = elem.getchildren()[0].getchildren()
                    if len(temp) > 0:
                        search_hit = temp[0]
                        search_hit = SearchHit( search_hit )
                        search_hit.set_spectrum_query( specQuery )
                        if getAll or (
                            search_hit.probability > threshold and not search_hit.decoy):
                            search_hits.append(  search_hit )
        return search_hits

    def get_threshold(self, myfdr):
        search_hits = self._parse_for_FDR( 0.0 )
        dec = [   s.probability for s in search_hits if s.decoy]
        nondec = [s.probability for s in search_hits if not s.decoy]
        hitratio = len( dec ) *1.0/ len( search_hits )
        #myfdr = 0.015
        fdr001 = -1
        breaks = 500
        for i in reversed(range( breaks )):
            n = len( [s for s in nondec if s >= i * 1.0 / breaks] )
            d = len( [s for s in dec if s >= i * 1.0 / breaks ] )
            fdr = 0
            if d + n > 0: fdr = 1.0-n*1.0/(d+n)  
            fdr_corrected = fdr / hitratio
            if fdr_corrected < myfdr: fdr001 = -1
            elif fdr001 == -1: fdr001 = i * 1.0 / breaks
        return fdr001
        #print i, n, d, fdr, fdr_corrected, fdr001

class Spectrum_Query:
    def __init__(self, elem):
        self.spectrum       = ( elem.get("spectrum") )
        self.charge         = int( elem.get("assumed_charge") )
        self.start_scan     = int( elem.get("start_scan") )
        self.end_scan       = int( elem.get("end_scan") )
        self.precursorMass  = float( elem.get("precursor_neutral_mass") )
        self.retTime        = -1
        if 'retention_time_sec' in elem.keys():
            self.retTime        = float( elem.get("retention_time_sec") )

    @property
    def scan_number(self):
        return self.start_scan 
        # myspec = self.spectrum.split('.')
        # assert self.charge == int(myspec[-1])
        # assert myspec[-2] == myspec[-3]
        # return int(myspec[-2])

class SearchHit:
    def __init__(self, search_hit):
        self.peptide = search_hit.get( "peptide" )
        self.protein = search_hit.get( "protein" )
        self.protein_descr = search_hit.get( "protein_descr" )
        self.matched_ions = search_hit.get( "num_matched_ions" )
        self.total_ions = search_hit.get( "tot_num_ions" )
        self.decoy = ( self.protein[:5] == 'DECOY' or 
                       self.protein[:7] == 'reverse' )
        self.massdiff = search_hit.get( "massdiff" )
        self.calcNeutralMass = float( search_hit.get( "calc_neutral_pep_mass") )
        self.modifications = []
        self.modified_peptide = self.peptide
        for child in search_hit.getchildren():
            if child.tag.find( "modification_info" ) != -1:
                self.modified_peptide = child.get( "modified_peptide" )
                for mod in child.getchildren():
                    self.modifications.append( [int(mod.get("position")),
                        float(mod.get("mass") ) ] )
            if child.get("name") == "dot": self.dot = float(child.get('value'))
            if child.get("name") == "delta": self.delta = float(child.get('value'))
            if child.get("name") == "dot_bias": self.db = float(child.get('value'))
            if child.get("name") == "fval": self.fval = float(child.get('value'))
            if child.get("name") == "pvalue": self.pvalue = float(child.get('value'))
            if child.get("name") == "precursor_mz_diff": 
                self.diff = float(child.get('value'))
            if child.tag.find( "analysis_result" ) != -1:
                for ana in child.getchildren():
                    if ana.tag.find( "peptideprophet_result" ) != -1:
                        self.probability = float(ana.get('probability'))
            #Mascot
            if child.get("name")=="ionscore"     :self.ionscore      = float(child.get('value'))
            if child.get("name")=="identityscore":self.identityscore = float(child.get('value'))
            if child.get("name")=="expect"       :self.expect        = float(child.get('value'))
        self.phospho = False
        self.phospho_modifications = []
        for mod in self.modifications:
            if self.peptide[ mod[0] -1  ] in ['S', 'T', 'Y']:
                self.phospho = True
                self.phospho_modifications.append( mod[0] -1 )
    def set_spectrum_query(self, spec):
        self.spectrum_query = spec

    @property
    def phospho_len(self): return len( self.phospho_modifications)

    @property
    def scan_number(self):
        return self.spectrum_query.scan_number

    @property
    def retention_time(self):
        res = self.spectrum_query.retTime
        if res == -1: raise ValueError("Retention time not in pepXML")
        return res

    def get_retention_time(self):
        return self.retention_time

    def get_assumed_charge(self):
        return self.spectrum_query.charge

    def get_neutral_precursor_mass(self):
        return self.spectrum_query.precursorMass

    def get_precursor_mass(self):
        return ( (self.spectrum_query.precursorMass +
            self.spectrum_query.charge ) / 
            self.spectrum_query.charge )

    def get_db_string(self, db, priority = -1):
        """Creates a string which will insert this hit into the db"""
        modstr = ''; modmass = 0
        for m in self.modifications:
            modstr += '%s, ' % (self.peptide[ m[0] -1 ])
            modmass += m[1]
        modstr = modstr[:-2]
        ch = self.get_assumed_charge()
        charged_mass = (self.calcNeutralMass + ch) / ch
        query = """insert into %s (sequence, 
                calc_neutral_mass, calc_charged_mass,
                modifications_string, modifications_mass, priority) 
                VALUES ('%s', %s, %s, '%s', %s, %s)""" % (db, self.peptide, 
                self.calcNeutralMass, charged_mass, modstr, modmass, priority)
        return query


