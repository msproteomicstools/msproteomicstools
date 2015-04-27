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
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""

class Peak:
    """Represents one peak of a spectrum."""
    def __init__(self, str=None, spectraST = False):
        self.spectraST = spectraST
        if str: self.parse_str( str )
    
    def initialize(self, peak, intensity, peak_annotation,statistics):
        self.peak            = float(peak)
        self.intensity       = float(intensity)
        self.peak_annotation = peak_annotation
        self.statistics      = statistics
        self._parse_statistics( self.statistics )
        self._parse_peak_annotation( self.peak_annotation , self.spectraST )
    
    def init_with_self( self, peak):
        import copy
        self.peak            = copy.deepcopy(peak.peak)
        self.intensity       = copy.deepcopy(peak.intensity)
        self.peak_annotation = copy.deepcopy(peak.peak_annotation)
        self.statistics      = copy.deepcopy(peak.statistics)
        self.spectraST         = copy.deepcopy(peak.spectraST)
        self._parse_statistics( self.statistics )
    
    def _parse_peak_annotation(self, peak_annotation, spectraST = False):
        '''By default, this below is the annotation used. It is not spectraST!!!'''
        #~ y2-NH3(24)
        #~ b-H2O-NH3(17)
        #~ a(17)
        #~ y(16)
        #~ y(16)_isotopicpeak
        #~ b-H2O(17)
        #~ b-NH3(17)
        #~ b(17)
        #~ b(17)_isotopicpeak
        #~ b+H2O(17)
        #~ y-H2O-NH3(17)
        #~ y-H2O(17)
        #~ y-NH3(17)
        #~ y(17)
        #~ y(17)_isotopicpeak
        #~ a-H2O(18)
        #~ a-NH3(18)
        #~ b-H2O-NH3(18)
        #~ a(18)
        #~ b-H2O(18)
        #~ b-NH3(18)
        #~ y-H2O-NH3(18)
        #~ b(18)
        #~ b(18)_isotopicpeak
        #~ y-H2O(18)
        #~ y-NH3(18)
        #~ b+H2O(18)
        #~ b-H2O-NH3(20)
        #~ y-H2O-NH3(20)
        #~ b+H2O(20)
        #~ a-H2O(21)
        #~ a-NH3(21)
        
        if spectraST : 
            self._parse_peak_annotation_spectraST(peak_annotation)
            return
        
        #raise "This is not spectraST!"
        
        self.frg_serie = ''
        self.frg_nr    = 0
        self.frg_z     = 1
        self.frg_is_isotope = False
        self.frg_loss = []
        self.is_frg_loss = False
        self.is_frg_gain = False
        self.frg_gain = []
        self.is_precursor = False
        self.is_unknown = peak_annotation[0] == '?'
        if self.is_unknown : return 
        
        #frg_serie
        self.frg_serie = peak_annotation[0]
        #frg_z
        if peak_annotation[1].isdigit() : self.frg_z = int(peak_annotation[1]) #only valid until charge state 9 (enough!)
        #isotope
        if len(peak_annotation)>13 and ( peak_annotation[-13:] == '_isotopicpeak' ) :
            self.frg_is_isotope = True
            peak_annotation = peak_annotation[:-13]
        
        #frg_nr & is_precursor
        spl_annotation = peak_annotation.split("(")
        if len(spl_annotation)>1 :
            self.frg_nr = int(spl_annotation[1][:-1])
        elif spl_annotation[0][:9] == 'precursor' : self.is_precursor = True
        
        #frg_loss
        spl_annotation = peak_annotation.split("-")
        if len(spl_annotation) > 1 :
            self.is_frg_loss = True
            for loss in spl_annotation[1:] :
                loss_clean = loss.split("(")[0]
                self.frg_loss.append ( loss_clean )
        
        #frg_gain
        spl_annotation = peak_annotation.split("+")
        if len(spl_annotation) > 1 :
            self.is_frg_gain = True
            for gain in spl_annotation[1:] :
                gain_clean = gain.split("(")[0]
                self.frg_gain.append ( gain_clean )
    
    
    def _parse_peak_annotation_spectraST(self, peak_annotation):
        '''This below is the spectraST annotation.'''
        #ITA/0.00,ITA/0.00
        #IIA/-0.00,IIA/-0.00,a1/-0.00
        #y3-18^2/0.06,y3-17^2/-0.44
        #y3^2/-0.99
        #b4-45^2/-0.00,b4-46^2/0.49
        #b2/-0.00,y2-46/-0.07,b4-44^2/0.48
        #b4-35^2/-0.95,y2-44/0.96
        #y2-35/0.00
        #b4-17^2/0.05,b4-18^2/0.55
        #y2-18/-0.00
        #y2-17/-0.00,y4-46^2/-0.48
        #y4-18^2/-0.50,y4-17^2/-0.99
        #y2/-0.00
        #b3-45/0.00,b5-34^2/0.02,b5-35^2/0.52
        #a5^2/-0.96
        #b3-34/0.00
        #b5^2/0.02,a3/-0.00
        #y5-44^2/0.99
        #b3-17/-0.00
        #y5-36^2/-0.04
        #y5-35^2/-0.45
        #b3/-0.00
        #[b3/0.0100]
        #b3i/-0.02
        #y5^2/0.05
        #b6-36^2/-0.50
        #b6-35^2/-0.02,b6-34^2/-0.52
        #a6^2/-0.50
        #b6^2/-0.52
        #y3-18/-0.01
        #y3-17/-0.95
        #y3/-0.00
        #b4-45/-0.01,b4-46/0.97,b7-34^2/0.01,b7-35^2/0.51,b7-36^2/1.00
        #b4-44/0.98,a7^2/-1.00
        #y4-17/0.04
        #b5-46/0.00
        #b5-45/-0.00
        #b5-46i/0.01
        #p-44^2/0.43,p-45^2/0.95
        #[y4/-0.0188]
        #y4/0.00,p-36^2/-0.50
        #p-35^2/-0.25,y4i/-0.25
        #p-35^2i/-0.53
        #p-35^2i/0.01,b5-35/-0.49
        #p-17^2/-0.02,p-18^2/0.47,a5/0.94
        #p^2/-0.00,b5-18/-0.51
        #p^2i/-0.00
        #b5-17/0.01
        #y5-44/-0.97
        #y5-35/0.02
        #y5-18/0.01,y5-17/-0.98
        #y5/0.01
        #b6/0.97
        #y6-18/-0.01,y6-17/-0.99
        #y6/-0.01
        #y6i/0.03
    
            
        self.frg_serie = ''
        self.frg_nr    = 0
        self.frg_z     = 1
        self.frg_is_isotope = False
        self.frg_loss = []
        self.is_frg_loss = False
        self.is_frg_gain = False
        self.frg_gain = []
        self.is_precursor = False
        self.is_unknown = peak_annotation[0] == '?'
        self.is_not_unique = False
        self.multiple_interpretation = False
        self.is_immonium = False
        self.mass_error  = 0.0
        
        if self.is_unknown : return 
    
        peak_annotations = peak_annotation.split(',')
        if len(peak_annotations) > 1 : self.multiple_interpretation    = True
    
        #So far, we pick only the first interpretation
        first_interp = str(peak_annotations[0])
        
        #is_not_unique
        if '[' in first_interp :
            self.is_not_unique = True
            first_interp = first_interp.strip("[]")
            
        if str(first_interp)[0] == 'I':
            self.is_immonium = True
            self.immonium_aa = first_interp[1]
            self.frg_serie   = 'immonium'
            if '/' in first_interp : self.mass_error = float(first_interp.split('/')[1])
            return
    
        
        #frg_serie
        self.frg_serie = str(first_interp)[0]
        
        #is_precursor
        if self.frg_serie == 'p' : self.is_precursor = True
        
        #frg_nr (it works until 99, enough!)
        if str(first_interp)[1].isdigit() : 
            if str(first_interp)[2].isdigit()    : 
                self.frg_nr = int(str(first_interp)[1:3])
            else                            : self.frg_nr = int(str(first_interp)[1])
        
        #frg_is_isotope
        if 'i' in first_interp : self.frg_is_isotope = True
        
        #frg_z (it works until charge state +9, enough!)
        if '^' in first_interp :
            self.frg_z = int(first_interp.split('^')[1][0])
        
        #frg_loss and frg_gain and mass precision
        annotation_ion = first_interp
        if '/' in annotation_ion : 
            self.mass_error = float(annotation_ion.split('/')[1])
            annotation_ion  = annotation_ion.split('/')[0]
            
        if '-' or '+' in annotation_ion :
            if '-' in annotation_ion : 
                frg_lg_txt = annotation_ion.split('-')[1]
                if 'i' in frg_lg_txt : frg_lg_txt = frg_lg_txt.split('i')[0]  #remove non desired symbols (isotope, charge state...)
                if '^' in frg_lg_txt : frg_lg_txt = frg_lg_txt.split('^')[0]  #remove non desired symbols (isotope, charge state...)
                self.frg_loss.append( float( frg_lg_txt ) )
                self.is_frg_loss = True
            if '+' in annotation_ion :
                frg_lg_txt = annotation_ion.split('+')[1]
                if 'i' in frg_lg_txt : frg_lg_txt = frg_lg_txt.split('i')[0]  #remove non desired symbols (isotope, charge state...)
                if '^' in frg_lg_txt : frg_lg_txt = frg_lg_txt.split('^')[0]  #remove non desired symbols (isotope, charge state...)
                self.frg_gain.append( float( frg_lg_txt ) )
                self.is_frg_gain = True
        
        #print "spectraST" , peak_annotations, first_interp, self.frg_serie , self.frg_nr , self.frg_z , self.frg_is_isotope , self.frg_loss , self.is_frg_loss , self.is_frg_gain , self.frg_gain , self.is_precursor , self.is_unknown , self.is_not_unique , self.multiple_interpretation , self.is_immonium , self.mass_error
        
    
    
    
    def _parse_statistics(self, stat):
        if stat.strip() == '':
            self.nr_replicates = self.nr_replicates_necessary = -1
            self.mzStdDev = self.intensityStdDevOverMean = -1
            return
        stripped = stat.strip()
        split = stripped.split(" ")
        s_lib_occurence = split[0].split("/")
        self.nr_replicates = int(s_lib_occurence[0])
        self.nr_replicates_necessary = int(s_lib_occurence[1])
        statistics = split[1].split("|")
        if len(statistics) > 0: self.mzStdDev = float(statistics[0])
        else: self.mzStdDev = -1
        if len(statistics) > 1: self.intensityStdDevOverMean = float(statistics[1])
        else: self.intensityStdDevOverMean = -1
    
    def parse_str(self, peak):
        speak = peak.split( "\t")
        peak = float( speak[0] )
        intensity = float( speak[1] )
        statistics = ''
        if len(speak) > 2: peak_annotation = speak[2]
        if len(speak) > 3: statistics = speak[3]
        self.initialize(peak, intensity, peak_annotation, statistics)
    
    def to_write_string(self):
        s = ''
        s += repr(self.peak) + '\t'
        s += repr(self.intensity) + '\t'
        s += (self.peak_annotation) + '\t'
        s += (self.statistics) + '\t'
        return s
