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

import string

# Isotope Modification
# 0 means no modification
# 1 means N15 (heavy nitrogen)
NOISOTOPEMODIFICATION = 0
N15_ISOTOPEMODIFICATION = 1

class Residues:
    """
    A class that contains information elements, amino acids and modifications.
    It stores mainly masse of these but also chemical formulas.

        The most commonly used properties are: 
            - Residues.average_elments : element weights 
            - Residues.monoisotopic_elments : element weights 
            - Residues.aa_codes : Three and One letter amino acid codes
            - Residues.aa_names : English names of the amino acids
            - Residues.aa_sum_formulas_text : Chemical formulas of all amino acids
            - Residues.aa_sum_formulas: Chemical formulas of all amino acids as hash
            - Residues.mass_xxx: monoisotopic masses of different compounds (NH3, H2O, CO, HPO4 etc)
            - Residues.average_data: average weight of amino acids
            - Residues.monoisotopic_data: monoisotopic weight of amino acids
            - Residues.monoisotopic_mod: monoisotopic modification data
            - Residues.mod_mapping: mapping of + notation to absolute weight notation (K[+8] to K[136])
            - Residues.Hydropathy: Hydropathy of amino acids (gravy scores)
            - TODO hydrophobicity of amino acids 
            - TODO basicity of amino acids 
            - TODO helicity of amino acids 
            - Residues.pI: pI of amino acids 
    """


    # http://www.sisweb.com/referenc/source/exactmaa.htm
    # http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    average_elements = {
        'H' : 1.007825   * 99.99/100 + 2.014102  * 0.015/100,
        'N' : 14.003074  * 99.63/100 + 15.000109 * 0.37/100,
        'O' : 15.994915  * 99.76/100 + 16.999131 * 0.038/100  + 17.999159 * 0.20/100,
        'C' : 12.000000  * 98.90/100 + 13.003355 * 1.10,
        'P' : 30.973763 
    }

    monoisotopic_elements = {
        'H'   : 1.007825032,
        'H2'  : 2.01410178,

        'C'   : 12.000000,
        'C13' : 13.00335484,

        'N'   : 14.003074005,
        'N15' : 15.000108898,

        'O'   : 15.994914620,
        'O17' : 16.999132,
        'O18' : 17.999161,

        'P'   : 30.973762,
        'S'   : 31.972071
    }

    aa_codes = {
        'A' : 'Ala',
        'R' : 'Arg',
        'N' : 'Asn',
        'D' : 'Asp',
        'C' : 'Cys',
        'E' : 'Glu',
        'Q' : 'Gln',
        'G' : 'Gly',
        'H' : 'His',
        'I' : 'Ile',
        'L' : 'Leu',
        'K' : 'Lys',
        'M' : 'Met',
        'F' : 'Phe',
        'P' : 'Pro',
        'S' : 'Ser',
        'T' : 'Thr',
        'W' : 'Trp',
        'Y' : 'Tyr',
        'V' : 'Val',
        'C[160]' : 'Cys+CAM',
        'M[147]' : 'Met+Ox',
    }
    aa_codes_rev = dict([(v,k) for k,v in aa_codes.iteritems()])

    aa_names = {
        'A': 'Alanine',
        'B': 'Aspartic Acid or Asparagine',
        'C': 'Cysteine',
        'c': 'Modified cysteine' ,
        'D': 'Aspartate',
        'E': 'Glutamate',
        'F': 'Phenylalanine',
        'G': 'Glycine',
        'H': 'Histidine',
        'I': 'Isoleucine',
        'K': 'Lysine',
        'k': 'Lys->Cys substitution and carbamidomethylation (903)',
        'L': 'Leucine',
        'M': 'Methionine',
        'm': 'Modified methionine' ,
        'N': 'Asparagine',
        'P': 'Proline',
        'Q': 'Glutamine',
        'R': 'Arginine',
        'S': 'Serine',
        'T': 'Threonine',
        'V': 'Valine',
        'W': 'Tryptophan',
        'X': 'Leucine/Isoleucine',
        'Y': 'Tyrosine',
        'Z': 'Glutamic acid'
        }

    aa_sum_formulas_text = {
        'A' : 'C3H5ON',
        'R' : 'C6H12ON4',
        'N' : 'C4H6O2N2',
        'D' : 'C4H5O3N',
        'C' : 'C3H5ONS',
        'E' : 'C5H7O3N',
        'Q' : 'C5H8O2N2',
        'G' : 'C2H3ON',
        'H' : 'C6H7ON3',
        'I' : 'C6H11ON',
        'L' : 'C6H11ON',
        'K' : 'C6H12ON2',
        'M' : 'C5H9ONS',
        'F' : 'C9H9ON',
        'P' : 'C5H7ON',
        'S' : 'C3H5O2N',
        'T' : 'C4H7O2N',
        'W' : 'C11H10ON2',
        'Y' : 'C9H9O2N',
        'V' : 'C5H9ON'
    }

    #from http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
    aa_sum_formulas = {
        'A' : { 'C' : 3,  'H' : 5  , 'O' : 1, 'N' : 1  },
        'R' : { 'C' : 6,  'H' : 12 , 'O' : 1, 'N' : 4  },
        'N' : { 'C' : 4,  'H' : 6  , 'O' : 2, 'N' : 2  },
        'D' : { 'C' : 4,  'H' : 5  , 'O' : 3, 'N' : 1  },
        'C' : { 'C' : 3,  'H' : 5  , 'O' : 1, 'N' : 1, 'S' : 1  },
        'E' : { 'C' : 5,  'H' : 7  , 'O' : 3, 'N' : 1  },
        'Q' : { 'C' : 5,  'H' : 8  , 'O' : 2, 'N' : 2  },
        'G' : { 'C' : 2,  'H' : 3  , 'O' : 1, 'N' : 1  },
        'H' : { 'C' : 6,  'H' : 7  , 'O' : 1, 'N' : 3  },
        'I' : { 'C' : 6,  'H' : 11 , 'O' : 1, 'N' : 1  },
        'L' : { 'C' : 6,  'H' : 11 , 'O' : 1, 'N' : 1  },
        'K' : { 'C' : 6,  'H' : 12 , 'O' : 1, 'N' : 2  },
        'M' : { 'C' : 5,  'H' : 9  , 'O' : 1, 'N' : 1, 'S' : 1   },
        'F' : { 'C' : 9,  'H' : 9  , 'O' : 1, 'N' : 1  },
        'P' : { 'C' : 5,  'H' : 7  , 'O' : 1, 'N' : 1  },
        'S' : { 'C' : 3,  'H' : 5  , 'O' : 2, 'N' : 1  },
        'T' : { 'C' : 4,  'H' : 7  , 'O' : 2, 'N' : 1  },
        'W' : { 'C' : 11, 'H' : 10 , 'O' : 1, 'N' : 2  },
        'Y' : { 'C' : 9,  'H' : 9  , 'O' : 2, 'N' : 1  },
        'V' : { 'C' : 5,  'H' : 9  , 'O' : 1, 'N' : 1  },
        'C[160]' : { 'C' : 3+2,  'H' : 5+3  , 'O' : 1+1, 'N' : 1+1, 'S' : 1  }, # + CAM = H(3) C(2) N O 
        'M[147]' : { 'C' : 5,  'H' : 9  , 'O' : 1+1, 'N' : 1, 'S' : 1   },
    }

    mass_H = monoisotopic_elements['H']
    mass_N = monoisotopic_elements['N']
    mass_O = monoisotopic_elements['O']
    mass_C = monoisotopic_elements['C']
    mass_S = monoisotopic_elements['S']
    mass_P = monoisotopic_elements['P']
    mass_NH2 = mass_N + 2*mass_H
    mass_NH3 = mass_N + 3*mass_H
    mass_CO =  mass_C +   mass_O
    mass_H2O = mass_O + 2*mass_H
    mass_OH =  mass_O +   mass_H
    mass_H3PO4 = mass_P + mass_O * 4 + mass_H * 3
    mass_H1PO4 = mass_P + mass_O * 4 + mass_H * 1
    mass_H1PO3 = mass_P + mass_O * 3 + mass_H * 1
    mass_CAM = 2* mass_C + 4*mass_H + mass_O + mass_N #CH2-CONH2

    mass_C13 = monoisotopic_elements['C13']
    mass_N15 = monoisotopic_elements['N15']
    mass_diffC13 = mass_C13 - mass_C
    mass_diffN15 = mass_N15 - mass_N

    #: average weight of amino acids
    average_data = { 
        'A': ('Alanine', 71.0788),
        'B': ('Aspartic Acid or Asparagine', 114.5962), 
        'C': ('Cysteine', 103.1448), 
        'c': ('Modified cysteine' , 160.1448), # Add 57
        'D': ('Aspartate', 115.0886),
        'E': ('Glutamate', 129.1155),
        'F': ('Phenylalanine', 147.1766),
        'G': ('Glycine', 57.0519),
        'H': ('Histidine', 137.1411),
        'I': ('Isoleucine', 113.1594),
        'K': ('Lysine', 128.1741),
        'k': ('Lys->Cys substitution and carbamidomethylation (903)', 128.09496 + 32.0219),
        'L': ('Leucine', 113.1594),
        'M': ('Methionine', 131.1986),
        'm': ('Modified methionine' , 147.1986), # add 16
        'N': ('Asparagine', 114.1038),
        'P': ('Proline', 97.1167),
        'Q': ('Glutamine', 128.1307),
        'R': ('Arginine', 156.1875),
        'S': ('Serine', 87.0782),
        'T': ('Threonine', 101.1051),
        'V': ('Valine', 99.1326),
        'W': ('Tryptophan', 186.2132),
        'X': ('Leucine/Isoleucine', 113.1594), # Can't distinguish leucine/isoleucine.
        'Y': ('Tyrosine', 163.176),
        'Z': ('Glutamic acid, or glutamine', 128),
        }
    
    #e.g. from http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
    # see also http://www.sbeams.org/svn/sbeams/trunk/sbeams/lib/perl/SBEAMS/Proteomics/AminoAcidModifications.pm
    monoisotopic_data = { 
        'A': ('Alanine', 71.03711),
        'B': ('Aspartic Acid or Asparagine', 114.04293), 
        'C': ('Cysteine', 103.00919), 
        'D': ('Aspartate', 115.02694),
        'E': ('Glutamate', 129.04259),
        'F': ('Phenylalanine', 147.06841),
        'G': ('Glycine', 57.02146),
        'H': ('Histidine', 137.05891),
        'I': ('Isoleucine', 113.08406),
        'K': ('Lysine', 128.09496),
        'L': ('Leucine', 113.08406),
        'M': ('Methionine', 131.04049),
        'N': ('Asparagine', 114.04293),
        'P': ('Proline', 97.05276),
        'Q': ('Glutamine', 128.05858),
        'R': ('Arginine', 156.10111),
        'S': ('Serine', 87.03203),
        'T': ('Threonine', 101.04768),
        'V': ('Valine', 99.06841),
        'W': ('Tryptophan', 186.07931),
        'X': ('Leucine/Isoleucine', 113.08406), 
        'Y': ('Tyrosine', 163.06333),
        'Z': ('Glutamic acid, or glutamine', 128.05858),
        }

    monoisotopic_mod = {
        'c': ('Modified cysteine', monoisotopic_data["C"][1] + mass_CAM - mass_H ), # CAM replaces H
        #'c': ('Modified cysteine' , 160.00919), # Add 57
        'C[160]': ('Modified cysteine', monoisotopic_data["C"][1] + mass_CAM - mass_H ), # CAM replaces H
        'k': ('Lys->Cys substitution and carbamidomethylation (903)', 128.09496 + 31.935685),
        'N[115]': ('Asparagine', monoisotopic_data["N"][1] - mass_N - mass_H + mass_O),
        #'m': ('Modified methionine', 147.04049), # add 16
        'm': ('Modified methionine', monoisotopic_data["M"][1] + mass_O), # oxygen
        'M[147]': ('Modified methionine', monoisotopic_data["M"][1] + mass_O), # oxygen
        # SILAC labels 
        'K[136]' : ('heavy Lysine',   monoisotopic_data["K"][1] + 8.014199), #UniMod:259
        'R[166]' : ('heavy Arginine', monoisotopic_data["R"][1] + 10.008269), #UniMod:267
        'R[162]' : ('heavy Arginine', monoisotopic_data["R"][1]  + 6*mass_diffC13), #UniMod:188
        'V[104]' : ('heavy Valine',   monoisotopic_data["V"][1] + 5*mass_diffC13), # no unimod
        'V[105]' : ('heavy Valine',   monoisotopic_data["V"][1] + 5*mass_diffC13 + mass_diffN15), # unimod 268
        # Pyro Unimod 27 and 28
        'E[111]': ('pyro Glutamate', 129.04259 - mass_O - 2*mass_H),
        'Q[111]': ('pyro Glutamine', 128.05858 - mass_O - 2*mass_H),
        # Unimod 385 # Pyro-carbamidomethyl as a delta from Carbamidomethyl-Cys
        'C[143]': ('Pyro-carbamidomethyl cysteine' , monoisotopic_data["C"][1] + mass_CAM - mass_H - 3*mass_H - mass_N),  
        # Phospho
        'S[166]': ('Phospho Serine', 87.03203 + mass_H1PO3),
        'S[167]': ('Phospho Serine', 87.03203 + mass_H1PO3),
        'T[181]': ('Phospho Threonine', 101.04768 + mass_H1PO3),
        'Y[243]': ('Phospho Tyrosine', 163.06333 + mass_H1PO3),
    }

    mod_mapping = {
        "K[+8]" : "K[136]",
        "R[+10]": "R[166]",
        "M[+16]": "M[147]",
        "N[-1]" : "N[115]",
        "C[+57]": "C[160]",
        "C[+40]": "C[160]",
        "R[+6]" : "R[162]",
        "V[+5]" : "V[104]",
        "V[+6]" : "R[105]",
        "S[+80]" : "S[167]",
        "T[+80]" : "T[181]",
        "Y[+80]" : "Y[243]",
    }

    monoisotopic_data.update(monoisotopic_mod)

    #C[169] 58       => ? 
    #C[152] 2        => ?
    #W[202] 23       => Oxidation?


    """
    http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
    GRAVY (Grand Average of Hydropathy) 
    The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values [9] of all the amino acids, divided by the number of residues in the sequence. 
    Amino acid scale: Hydropathicity.

    Author(s): Kyte J., Doolittle R.F.
    Reference: J. Mol. Biol. 157:105-132(1982).

    Amino acid scale values:
    """


    # source http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
    Hydropathy = {
    'Ala':  1.800,  
    'Arg': -4.500,  
    'Asn': -3.500,  
    'Asp': -3.500,  
    'Cys':  2.500,  
    'Gln': -3.500,  
    'Glu': -3.500,  
    'Gly': -0.400,  
    'His': -3.200,  
    'Ile':  4.500,  
    'Leu':  3.800,  
    'Lys': -3.900,  
    'Met':  1.900,  
    'Phe':  2.800,  
    'Pro': -1.600,  
    'Ser': -0.800,  
    'Thr': -0.700,  
    'Trp': -0.900,  
    'Tyr': -1.300,  
    'Val':  4.200,  
    }
    Hydropathy_aa = dict([ (aa_codes_rev[k],v) for k,v in Hydropathy.iteritems()])

    # http://www.anaspec.com/html/pK_n_pl_Values_of_AminoAcids.html
    pI = {
        'G': 6.0,
        'A': 6.0,
        'V': 6.0,
        'L': 6.0,
        'X': 6.0, #L or I
        'I': 6.0,
        'F': 5.5,
        'P': 6.3,
        'S': 5.7,
        'T': 5.6,
        'Y': 5.7,
        'C': 5.0,
        'M': 5.7,
        'N': 5.4,
        'B': 4.1, #avg N and D
        'Q': 5.7,
        'Z': 4.5, #avg Q,E
        'W': 5.9,
        'D': 2.8,
        'E': 3.2,
        'K': 9.7,
        'R': 10.8,
        'H': 7.6
    }

    def __init__(self, type="mono"):
        #add the phosphorylations
        self.monoisotopic_data[ 's' ] = ('Phospho-S', 
        self.monoisotopic_data[ 'S' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 't' ] = ('Phospho-T', 
        self.monoisotopic_data[ 'T' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 'y' ] = ('Phospho-Y',
        self.monoisotopic_data[ 'Y' ][1] + self.mass_H1PO3)

        self.average_data[ 's' ] = ('Phospho-S', 
        self.average_data[ 'S' ][1] + self.mass_H1PO3)
        self.average_data[ 't' ] = ('Phospho-T', 
        self.average_data[ 'T' ][1] + self.mass_H1PO3)
        self.average_data[ 'y' ] = ('Phospho-Y',
        self.average_data[ 'Y' ][1] + self.mass_H1PO3)

        if not type:
            self.residues = self.average_data
        elif type.startswith("mono"):
            self.residues = self.monoisotopic_data
        elif type.startswith("av"):
            self.residues = self.average_data
        else:
            raise ValueError("Type of residue must be one of: mono[isotopic], av[erage] (characters within [] are optional.")
        keys = self.residues.keys()
        self.res_pairs = [ string.join((r, s), '') for r in keys for s in keys ]

    def recalculate_monisotopic_data(self):

        self.monoisotopic_data = {}
        for abbrev, formula in self.aa_sum_formulas.iteritems(): 
            mysum = 0.0
            for key, value in formula.iteritems():
                mysum += self.monoisotopic_elements[ key ] * value
            self.monoisotopic_data[ abbrev ] = ( self.aa_codes[abbrev] , mysum )

        #
        self.monoisotopic_data['c'] = self.monoisotopic_data['C'] + self.mass_CAM - self.mass_H
        self.monoisotopic_data['c'] = ( 'Modified cystein', 
             self.monoisotopic_data['C'][1] + self.mass_CAM - self.mass_H)
        self.monoisotopic_data['k'] = ( 'Lys->Cys substitution and carbamidomethylation (903)',
            self.monoisotopic_data['K'][1] + 31.935685)
        self.monoisotopic_data['m'] = ( 'Modified methionine', 
             self.monoisotopic_data['M'][1] + self.mass_O)
        self.monoisotopic_data[ 's' ] = ('Phospho-S', 
            self.monoisotopic_data[ 'S' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 't' ] = ('Phospho-T', 
            self.monoisotopic_data[ 'T' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 'y' ] = ('Phospho-Y',
            self.monoisotopic_data[ 'Y' ][1] + self.mass_H1PO3)
        self.residues = self.monoisotopic_data

    def recalculate_monisotopic_data_for_N15(self):

        self.monoisotopic_data = {}
        for abbrev, formula in self.aa_sum_formulas.iteritems(): 
            mysum = 0.0
            for key, value in formula.iteritems():
                #replace N with N15
                if key == 'N': key = 'N15'
                mysum += self.monoisotopic_elements[ key ] * value
            self.monoisotopic_data[ abbrev ] = ( self.aa_codes[abbrev] , mysum )

        #IMPORTANT: CAM is added afterwards and is NOT heavy
        #
        self.monoisotopic_data['C[160]'] = ( 'Modified cystein', 
             self.monoisotopic_data['C'][1] + self.mass_CAM - self.mass_H)
        self.monoisotopic_data['N[115]'] = ( 'Modified asparagine', 
             self.monoisotopic_data['N'][1] - self.mass_N15 - self.mass_H + self.mass_O) 
        self.monoisotopic_data['M[147]'] = ( 'Modified methionine', 
             self.monoisotopic_data['M'][1] + self.mass_O)
        #
        self.monoisotopic_data['c'] = ( 'Modified cystein', 
             self.monoisotopic_data['C'][1] + self.mass_CAM - self.mass_H)
        self.monoisotopic_data['k'] = ( 'Lys->Cys substitution and carbamidomethylation (903)',
            self.monoisotopic_data['K'][1] + 31.935685)
        self.monoisotopic_data['m'] = ( 'Modified methionine', 
             self.monoisotopic_data['M'][1] + self.mass_O)
        self.monoisotopic_data[ 's' ] = ('Phospho-S', 
            self.monoisotopic_data[ 'S' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 't' ] = ('Phospho-T', 
            self.monoisotopic_data[ 'T' ][1] + self.mass_H1PO3)
        self.monoisotopic_data[ 'y' ] = ('Phospho-Y',
            self.monoisotopic_data[ 'Y' ][1] + self.mass_H1PO3)
        self.residues = self.monoisotopic_data

