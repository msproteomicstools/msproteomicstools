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

import sys

class Elements:

    def __init__(self):
        self.list = []
        self._initElements()

    def _initElements(self):

        elH = Element('H',[1.007825032,2.014101778],[0.99984426,0.00015574])
        elC = Element('C',[12.000000000,13.00335484],[0.988922,0.011078])
        elN = Element('N',[14.00307401,15.0001089],[0.996337,0.003663])
        elO = Element('O',[15.99491462,16.9991315,17.9991604],[0.997628,0.000372,0.002000])
        elS = Element('S',[31.97207069,32.9714585,33.96786683,35.96708088],[0.95018,0.0075,0.04215,0.00017])
        elP = Element('P',[30.97376151],[1.00000])


        self.list.append(elC)
        self.list.append(elH)
        self.list.append(elN)
        self.list.append(elO)
        self.list.append(elS)
        self.list.append(elP)
        
        
    def addElement(self,symbol,isotMass,isotAbundance):
        newEl = Element(symbol,isotMass)
        self.list.append(newEl)
        

class Element:
    
    def __init__(self,symbol,isotMass,isotAbundance):
        self.symbol = ''
        self.isotMass = []
        self.isotAbundance = []
        self.symbol = symbol
        self.isotMass = isotMass
        self.isotAbundance = isotAbundance
        
def test():
    isots = Elements()
    
    for el in isots.list:
        print el.symbol
        totalAb = 0
        for k,m in zip(el.isotMass,el.isotAbundance):
            print "Mass: %s, Abundance: %s" % (k,m)
            totalAb += m
        print "Total abundance : " , totalAb
            
    #Monoisotopic masses
    for el in isots.list:
        print el.symbol , el.isotMass[0]
        
    
if __name__=="__main__":
    test()
    sys.exit(2)