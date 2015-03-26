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
try:
    from elements import Elements
    from elements import Formulas
except ImportError:
    # Python 3
    from .elements import Elements
    from .elements import Formulas


class Aminoacides:
    
    def __init__(self):
        self.list = []
        self.initAminoacides()


    def initAminoacides(self):
        
        ala = Aminoacid('Alanine','A','Ala',{'C':3,'H':5,'N':1,'O':1})
        arg = Aminoacid('Arginine','R','Arg',{'C':6,'H':12,'N':4,'O':1})
        asn = Aminoacid('Asparagine','N','Asn',{'C':4,'H':6,'N':2,'O':2})
        asp = Aminoacid('Aspartic acid','D','Asp',{'C':4,'H':5,'N':1,'O':3})
        cys = Aminoacid('Cysteine','C','Cys',{'C':3,'H':5,'N':1,'O':1,'S':1})
        glu = Aminoacid('Glutamic acid','E','Glu',{'C':5,'H':7,'N':1,'O':3})
        gln = Aminoacid('Glutamine','Q','Gln',{'C':5,'H':8,'N':2,'O':2})
        gly = Aminoacid('Glycine','G','Gly',{'C':2,'H':3,'N':1,'O':1})
        his = Aminoacid('Histidine','H','His',{'C':6,'H':7,'N':3,'O':1})
        ile = Aminoacid('Isoleucine','I','Ile',{'C':6,'H':11,'N':1,'O':1})
        leu = Aminoacid('Leucine','L','Leu',{'C':6,'H':11,'N':1,'O':1})
        lys = Aminoacid('Lysine','K','Lys',{'C':6,'H':12,'N':2,'O':1})
        met = Aminoacid('Methionine','M','Met',{'C':5,'H':9,'N':1,'O':1,'S':1})
        phe = Aminoacid('Phenylalanine','F','Phe',{'C':9,'H':9,'N':1,'O':1})
        pro = Aminoacid('Proline','P','Pro',{'C':5,'H':7,'N':1,'O':1})
        ser = Aminoacid('Serine','S','Ser',{'C':3,'H':5,'N':1,'O':2})
        thr = Aminoacid('Threonine','T','Thr',{'C':4,'H':7,'N':1,'O':2})
        trp = Aminoacid('Tryptophan','W','Trp',{'C':11,'H':10,'N':2,'O':1})
        tyr = Aminoacid('Tyrosine','Y','tyr',{'C':9,'H':9,'N':1,'O':2})
        val = Aminoacid('Valine','V','Val',{'C':5,'H':9,'N':1,'O':1})
        
        self.addAminoacid(ala)
        self.addAminoacid(arg)
        self.addAminoacid(asn)
        self.addAminoacid(asp)
        self.addAminoacid(cys)
        self.addAminoacid(glu)
        self.addAminoacid(gln)
        self.addAminoacid(gly)
        self.addAminoacid(his)
        self.addAminoacid(ile)
        self.addAminoacid(leu)
        self.addAminoacid(lys)
        self.addAminoacid(met)
        self.addAminoacid(phe)
        self.addAminoacid(pro)
        self.addAminoacid(ser)
        self.addAminoacid(thr)
        self.addAminoacid(trp)
        self.addAminoacid(tyr)
        self.addAminoacid(val)
        
    def addAminoacid(self,aminoacid):
        self.list.append(aminoacid)

    def getAminoacid(self,code):
        for aa in self.list:
            if aa.code == code: return aa
        raise Exception('Element %s is not implemented' % code)
            
class Aminoacid:
    """
    Class to hold information about a single Amino Acid (AA)
    """
    
    def __init__(self,name, code, code3,composition):
        #: Full name of the AA
        self.name = name
        #: One letter code
        self.code = code
        #: Three letter code
        self.code3 = code3
        #: Elemental composition
        self.composition = {}
        self.composition = composition
        #: Library of elements
        self.elementsLib = Elements()
        self.deltaMass = Formulas.mass(self.composition, elementsLib = self.elementsLib) #self.calDeltaMass()
    
    
if __name__ == "__main__":
    
    myAAs = Aminoacides()
    
#    for aa in myAAs.list:
#        for elem,num in aa.composition.iteritems():
#            print("%s: %s" % (elem,num))
    
    print("1-letter code    3-letter code    Chemical formula    Monoisotopic")
    print("-"*65)
    for aa in myAAs.list:
        print("%-10s  %-15s  %-20s %s" % (aa.code, aa.code3, Formulas.compositionString(aa.composition), aa.deltaMass))
          
##1-letter code 3-letter code Chemical formula Monoisotopic
##-----------------------------------------------------------------
##A           Ala              H5C3ON               71.03711379
##R           Arg              H12C6ON4             156.101111044
##N           Asn              H6C4O2N2             114.042927452
##D           Asp              H5C4O3N              115.02694303
##C           Cys              H5C3SON              103.00918448
##E           Glu              H7C5O3N              129.042593094
##Q           Gln              H8C5O2N2             128.058577516
##G           Gly              H3C2ON               57.021463726
##H           His              H7C6ON3              137.058911874
##I           Ile              H11C6ON              113.084063982
##L           Leu              H11C6ON              113.084063982
##K           Lys              H12C6ON2             128.094963024
##M           Met              H9C5SON              131.040484608
##F           Phe              H9C9ON               147.068413918
##P           Pro              H7C5ON               97.052763854
##S           Ser              H5C3O2N              87.03202841
##T           Thr              H7C4O2N              101.047678474
##W           Trp              H10C11ON2            186.07931296
##Y           tyr              H9C9O2N              163.063328538
##V           Val              H9C5ON               99.068413918  
