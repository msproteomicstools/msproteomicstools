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

from __future__ import print_function
import sys
import re

class Formulas :
    
    #This mass tolerance is used only for matching isotope masses! 
    #Since isotopes are typically reported as rounded masses (i.e. C13), this value should be quite high
    massTolerance_forMatch = 0.3   
    
    
    H2O         = { 'H' : 2 , 'O' : 1}
    NH2         = { 'H' : 2 , 'N' : 1}
    NH3         = { 'H' : 3 , 'N' : 1}
    CO          = { 'C' : 1 , 'O' : 1}
    CO2         = { 'C' : 1 , 'O' : 2}
    OH          = { 'H' : 1 , 'O' : 1}
    H3PO4       = { 'H' : 3 , 'P' : 1, 'O' : 4 }
    HPO4        = { 'H' : 1 , 'P' : 1, 'O' : 4 }
    HPO3        = { 'H' : 1 , 'P' : 1, 'O' : 3 }
    CH2_CONH2   = { 'H' : 4 , 'N' : 1, 'C' : 2, 'O' : 1 }  #CH2-CONH2 (CAM)
    
    def __init__(self) :
        pass
        
    @staticmethod
    def mass(formula, elementsLib = None):
        elements = elementsLib
        if elementsLib == None : elements = Elements()
        deltaMass=0
        
        for el , numatoms in formula.items() : 
            #Search it in the elements list
            found = False
            for el2 in elements.list :
                if el2.symbol == el : deltaMass += el2.isotMass[0] * numatoms ; found = True
            if not found : #it *might* be an isotope of an element --> figure it out
                element_match = re.findall('[A-Za-z]+' , el )
                isotope_match = (re.findall('[\d\.]+', el))
                if not len(isotope_match) == 1 or not len(element_match) == 1 : 
                    #Throw an Exception. Note : Should we create an Exception handler object??
                    print("The following composition has not been recognized : " , formula)
                    print("Please, review the compositions of your modifications as a probable cause of this error.")
                    print("Error caused by : " , element_match , isotope_match)
                    sys.exit(5)
                
                found2 = False
                for el in elements.list :
                    if element_match[0] == el.symbol : 
                        #find the closest isotope match for this element
                        for isotopeMass in el.isotMass : 
                            if abs(isotopeMass - float(isotope_match[0])) < Formulas.massTolerance_forMatch : 
                                deltaMass += isotopeMass * numatoms ; 
                                found2 = True
                if not found2 : 
                    #Throw an Exception. Note : Should we create an Exception handler object??
                    print("The following composition has not been recognized : " , formula, deltaMass)
                    print("Please, review the compositions of your modifications as a probable cause of this error.")
                    print("Error caused by : " , element_match[0] , isotope_match[0])
                    sys.exit(5)                        
                
        return deltaMass

    @staticmethod
    def add2components(formula1,formula2) :
        formula = formula1.copy()
        for el, value in formula2.items() :
            if el in formula    : formula[el] = formula[el] + formula2[el]
            else                : formula[el] = formula2[el]
        return formula
    
    @staticmethod
    def substract2components(formula1,formula2) :
        formula = formula1.copy()
        for el,value in formula2.items() :
            if el in formula : formula[el] = formula[el] - formula2[el]
            else             : formula[el] = -formula2[el]
        return formula

    @staticmethod
    def compositionString(formula):
        compString=""
        for elem,num in sorted(formula.items()):
            if num>1:
                compString += "%s%s" % (elem,num)
            else:
                compString += "%s" % elem
        return compString

class Elements:

    def __init__(self):
        self.list = []
        self._initElements()
        self.element = {}
        for el in self.list :
            self.element[el.symbol] = el
            
    def _initElements(self):

        self.addElement('H',[1.007825032,2.014101778],[0.99984426,0.00015574])
        self.addElement('C',[12.000000000,13.00335484],[0.988922,0.011078])
        self.addElement('N',[14.00307401,15.0001089],[0.996337,0.003663])
        self.addElement('O',[15.99491462,16.9991315,17.9991604],[0.997628,0.000372,0.002000])
        self.addElement('S',[31.97207069,32.9714585,33.96786683,35.96708088],[0.95018,0.0075,0.04215,0.00017])
        self.addElement('P',[30.97376151],[1.00000])
        
    def addElement(self,symbol,isotMass,isotAbundance):
        # TODO this is broken, Elements signature does not match
        
        #Check that the lists of masses and natural abundances have the same length
        if len(isotMass) != len(isotAbundance) :
            #Throw an exception
            print("Error : the isotopic masses and the natural abundance vector sizes don't match!")
            print("Element : " , symbol)
            print("isotopic masses : " , isotMass)
            print("natural abundances : " , isotAbundance)
            sys.exit(5)

        #Check that the sum of all the natural abundances is close enough to 1 (over 1 is not good, a bit below 1 might be acceptable)
        sumAbundances = sum(isotAbundance)
        
        if sumAbundances > 1 : 
            print("Error : the sum of the abundances is over 1!")
            print("Element : " , symbol)
            print("isotopic masses : " , isotMass)
            print("natural abundances : %s , that makes : %s" % (isotAbundance, sumAbundances))
            sys.exit(5)
        
        if sumAbundances < 0.97 :
            print("Error : the sum of the abundances is too low! It should be closer to 1.")
            print("Element : " , symbol)
            print("isotopic masses : " , isotMass)
            print("natural abundances : %s , that makes : %s" % (isotAbundance, sumAbundances))
            sys.exit(5)
            
        newEl = Element(symbol,isotMass, isotAbundance)
        self.list.append(newEl)

    def getElement(self,symbol):
        for el in self.list:
            if el.symbol == symbol: return el
        raise Exception('Element %s is not implemented' % symbol)
        

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
        print (el.symbol)
        totalAb = 0
        for k,m in zip(el.isotMass,el.isotAbundance):
            print("Mass: %s, Abundance: %s" % (k,m))
            totalAb += m
        print("Total abundance : " , totalAb)
            
    #Monoisotopic masses
    for el in isots.list:
        print(el.symbol , el.isotMass[0])
    
    print("Adding H2O and H3PO4 results in :", Formulas.add2components(Formulas.H2O, Formulas.H3PO4))
    print("Substracting H2O to H3PO4 results in :", Formulas.substract2components(Formulas.H3PO4, Formulas.H2O))
    print("Calculating the mass for the formula 13C6_18O2 : " , Formulas.mass({'C13' : 6 , '18O' : 2}))
    #print (isots.element)
    
if __name__=="__main__":
    test()
    sys.exit(2)
