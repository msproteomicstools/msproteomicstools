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
import re

"""
Doc :
    Contains some classes that can be useful in proteomics.
    Organism
    Protein - can be trypsinized
    Peptide - can be fragmented
    Fragment
"""

class Organism:
    def __init__(self):
        pass
    def add_proteins(self, proteins):
        self.seq_dic = {}
        for p in proteins:
            for pep in p.peptides:
                if self.seq_dic.has_key( pep.sequence): 
                    self.seq_dic[ pep.sequence ].append( pep )
                else: self.seq_dic[ pep.sequence ] = [pep]
    def find_peptides_from_same_protein(self, allseqs):
        queryseqs = []
        for seq in allseqs:
            try:
                prot = [pep.protein for pep in self.seq_dic[ seq ] ]
                for p in prot:
                    queryseqs.extend( [pep.sequence for pep in p.peptides if
                        pep.genome_occurence == 1 and not pep.sequence == seq] )
            except KeyError: pass
        return queryseqs
    def find_proteins_of_peptides(self, sequences):
        result = {}
        for pep in sequences:
            result[ org.seq_dic[ pep ][0].protein.id ] = -1
        return result

class Protein:
    def __init__(self):
        pass
    def initialize_dbrows(self, t, rows):
        self.id = t.row(rows[0], 'protein_key')
        self.sequence_key =  t.row(rows[0], 'sequence_key')
        self.experiment_key =  t.row(rows[0], 'experiment_key')
        self.peptides = []
        for row in rows:
            p = Peptide()
            p.initialize_dbrow( t, row )
            p.add_protein( self )
            self.peptides.append( p )

    @staticmethod
    def trypsin(aminoacids):
        sub = ''
        while aminoacids:
            k, r = aminoacids.find('K'), aminoacids.find('R')
            if k >= 0 and r >= 0 : cut = min(k, r)+1 
            else: cut =  max(k, r)+1
            # cut = min(k, r)+1 if k >= 0 and r >= 0 else max(k, r)+1
            if max(k,r) == -1: cut = len(aminoacids)
            sub += aminoacids[:cut]
            aminoacids = aminoacids[cut:]
            if not aminoacids or aminoacids[0] != 'P':
                yield sub 
                sub = ''

class Peptide:

    def __init__(self):
        self._raw_sequence  = None

    def __len__(self):
        return len( self.get_raw_sequence() )

    def init_with_self(self, other):
        self.set_sequence(other.sequence)
        self.ssr_calc  = other.ssr_calc
        self.id  = other.id

    def initialize_dbrow(self, t, row):
        self.peptide_key =  t.row(row, 'peptide_key')
        self._sequence =  t.row(row, 'sequence')
        self.molecular_weight =  t.row(row, 'molecular_weight')
        #try: 
        self.genome_occurence =  t.row(row, 'genome_occurence')
        #except KeyError: pass

    def add_protein(self, protein):
        self.protein = protein

    def set_sequence( self, sequence, format = 'bracket' ):
        self._raw_sequence = None
        if format == 'bracket': self._sequence = sequence
        elif format == 'SEQUEST':
            #re.sub( '@', '\[147\]', \
            #re.sub( '#', '\[160\]', \
            #re.sub( '#', 'c\[31\]', \
            self._sequence = \
            re.sub( 'S\*', 'S[167]', 
            re.sub( 'T\*', 'T[181]', 
            re.sub( 'Y\*', 'Y[243]', 
                    sequence ))) #)))
        else: raise ValueError('cannot set other format than bracket and SEQUEST')

    @property
    def sequence(self):
        return self.get_raw_sequence()

    def get_raw_sequence(self):
        if self._raw_sequence is None:
            self._raw_sequence = re.sub( '[^A-Z]', '', self._sequence )
        return self._raw_sequence
    #sequence = property(  get_raw_sequence, set_sequence )

    def get_modified_sequence(self, format = 'bracket' ):
        if format == 'bracket': return self._sequence 
        elif format == 'SEQUEST': 
            return \
            re.sub( '\[147\]', '@', 
            re.sub( '\[160\]', '#', 
            re.sub( 'c\[31\]', '#', 
            re.sub( '\[167\]', '*', 
            re.sub( '\[181\]', '*', 
            re.sub( '\[243\]', '*', 
                    self._sequence ))) )))
        else: raise ValueError('format %s not supported' % format)

    def has_phospho(self):
        tmp = self.get_modified_sequence( 'SEQUEST' )
        if tmp.find( '*' ) != -1: return True
        return False

    def modify_cysteins(self):
        #Alkylate Cysteins with CAM
        if self._sequence.find( 'C[') != -1: 
            raise AssertionError('already modified cysteins here!')
        self._sequence = self._sequence.replace( 'C', 'C[160]')

    def oxidize_methionines(self):
        # Oxidize Methionines
        if self._sequence.find( 'M[') != -1: 
            raise AssertionError('already oxidized methionines here!')
        self._sequence = self._sequence.replace( 'M', 'M[147]')

    def get_phospho_position(self):
        modseq = self.get_modified_sequence('SEQUEST')
        pos = []
        for i,aa in enumerate(modseq):
            #subtract the number of stars we already passed
            #subtract 1 because we start with lenght 0
            if aa == '*': pos.append( i - 1 - len( pos ) )
        return pos

    def get_maximal_charge(self):
        """ Count the number of amino acids that can hold a charge: R (Arg), H
        (His) or K (Lys) and add 1 for the N-terminus"""
        return self.sequence.count('R') + \
               self.sequence.count('H') + \
               self.sequence.count('K') + 1

    def create_fragmentation_pattern(self, R, bions=True, yions=True,
             aions=False, aMinusNH3=False, 
             bMinusH2O=False, bMinusNH3=False, bPlusH2O=False,
             yMinusH2O=False, yMinusNH3=False, cions=False, xions=False,
             zions=False, MMinusH2O=False, MMinusNH3=False ):
        seq = self.get_modified_sequence('bracket')
        self.mass = 0
        fragment_series = []
        #each fragment mass is an element of this species:
        #                    
        #                   R
        #                   | 
        #                HN-C-C=O
        #                   |
        #                   H
        #
        # 
        # The b series (from N terminus) looks like this (triple bond to O)
        #                    
        #                   R
        #                   | 
        #               H2N-C-C\=O
        #                   |
        #                   H
        #
        #
        # The y series (from C terminus) looks like this
        #                    
        #                   R
        #                   | 
        #              +H3N-C-COOH
        #                   |
        #                   H
        #
        #
        # see also http://www.matrixscience.com/help/fragmentation_help.html
        # note that the b and y series only go up to y[n-1] and b[n-1] since
        # the last ion of the series would be the parent ion (y) or the parent
        # ion with a loss of water (b).
        for q in re.finditer( '([A-Z]\[\d*\]|[A-Z])', seq):
            element = q.group(0)
            res_mass = R.residues[element][1]
            self.mass += res_mass
            fragment_series.append( self.mass )

        #we add a Water H2O to the complete peptide and to the charged peptide
        #as many protons as there are charges. 
        #To the b series we add one H because we want to start with NH2 and not NH
        #To the y series we add one OH and two protons since we want to end with
        #COOH and not CO and we start with NH3+ and not HN
        self.b_series = [b + R.mass_H for b in fragment_series[:-1]] 
        self.y_series = [self.mass - y + 2*R.mass_H + R.mass_OH for y in fragment_series[:-1]]
        ### print fragment_series
        # 
        self.allseries = []
        if bions: self.allseries.extend(self.b_series)
        if yions: self.allseries.extend(self.y_series)
        #these are all the more exotic ions
        if aions: self.a_series = [f - R.mass_CO for f in self.b_series]; self.allseries.extend(self.a_series)
        if cions: self.c_series = [f + R.mass_NH3 for f in self.b_series]; self.allseries.extend(self.c_series)
        if xions: self.x_series = [f + R.mass_CO - 2*R.mass_H for f in self.y_series]; self.allseries.extend(self.x_series)
        if zions: self.z_series = [f - R.mass_NH3 for f in self.y_series]; self.allseries.extend(self.z_series)
        if aMinusNH3: self.a_minus_NH3 = [f - R.mass_CO - R.mass_NH3 for f in self.b_series]; self.allseries.extend(self.a_minus_NH3)
        if bMinusH2O: self.b_minus_H2O = [b - R.mass_H2O for b in self.b_series]; self.allseries.extend(self.b_minus_H2O)
        if bMinusNH3: self.b_minus_NH3 = [b - R.mass_NH3 for b in self.b_series]; self.allseries.extend(self.b_minus_NH3)
        if bPlusH2O:  self.b_plus_H2O  = [b + R.mass_H2O for b in self.b_series]; self.allseries.extend(self.b_plus_H2O)
        if yMinusH2O: self.y_minus_H2O = [y - R.mass_H2O for y in self.y_series]; self.allseries.extend(self.y_minus_H2O)
        if yMinusNH3: self.y_minus_NH3 = [y - R.mass_NH3 for y in self.y_series]; self.allseries.extend(self.y_minus_NH3)
        if MMinusH2O: self.waterloss = self.mass; self.allseries.append(self.waterloss)
        if MMinusNH3: self.nh3loss =   self.mass + R.mass_H2O - R.mass_NH3; self.allseries.append(self.nh3loss)

        #
        self.mass += R.mass_OH + R.mass_H 

        self.molecular_weight = self.mass
        if self.charge > 0:
          self.charged_mass = (self.mass + self.charge * R.mass_H) / self.charge
        del self.mass

    def missed_cleavages(self):
        count = 0
        seq = self.get_raw_sequence()
        for i in range(len(seq)-1):
            if (seq[i] == 'K' or seq[i] == 'R') and not \
                seq[i+1] == 'P':
                count += 1
        return count

    def _get_modified_fragments(self):
        seq = self.get_modified_sequence('bracket')
        for q in re.finditer( '([A-Z]\[\d*\]|[A-Z])', seq):
            yield q.group(0)

    def get_fragment_objects(self, fragments, series, charge, R, q3_low, q3_high):

            ch = charge
            scount = 0
            for pred in fragments:
                scount += 1
                q3 = ( pred + (ch -1)*R.mass_H)/ch
                if q3 < q3_low or q3 > q3_high: continue
                yield Fragment(q3, series + str(scount), ch )

    def get_GRAVY(self):
        """
        http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
        GRAVY (Grand Average of Hydropathy) 
        The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values [9] of all the amino acids, divided by the number of residues in the sequence. 
        Amino acid scale: Hydropathicity.

        Author(s): Kyte J., Doolittle R.F.
        Reference: J. Mol. Biol. 157:105-132(1982).
        """

        import Residues
        return sum([ Residues.Residues.Hydropathy_aa[aa] for aa in self.sequence]) / len(self.sequence)

class Fragment():

    def __init__(self, q3=None, annotation=None, charge=None): 
        self.q3             = q3
        self.annotation     = annotation
        self.charge         = charge

    def __repr__(self):
        return "Fragment: %s (%s) %s+" %(round(self.q3, 2), self.annotation, self.charge)

    @property
    def pQ3(self):
        return round(self.q3, 2)


def test_fragmentation():
    import Residues
    R = Residues.Residues('mono')
    peptide = Peptide()
    peptide.charge = 1
    peptide.set_sequence('AAAAEIAVK')
    peptide.create_fragmentation_pattern(R)
    b = peptide.b_series
    y = peptide.y_series
    #
    #values from spectral library
    assert abs(y[-2] - 246.4 + 0.22) < 1e-2
    assert abs(b[4-1] - 285.3 + 0.14) < 1e-2
    assert abs(y[-3] - 317.1 - 0.12) < 1e-2
    assert abs(y[-6] - 630.4 + 0.02) < 1e-2
    assert abs(b[8-1] - 697.2 - 0.19) < 1e-2

    peptide = Peptide()
    peptide.charge = 2
    peptide.set_sequence('AAAAEIAVK')
    peptide.create_fragmentation_pattern(R)
    b = peptide.b_series
    y = peptide.y_series
    ch = peptide.charge
    y2 =  [ ( pred + (ch -1)*R.mass_H)/ch for pred in y]
    b2 =  [ ( pred + (ch -1)*R.mass_H)/ch for pred in b]
    #
    #values from spectral library
    assert abs(y2[-3] - 159.6 + 0.49) < 1e-2
    assert abs(b[3-1] - 214.0 - 0.12) < 1e-2
    assert abs(y[-2] - 246.2 + 0.02) < 1e-2
    assert abs(y2[-6] - 314.8 - 0.89) < 1e-2
    assert abs(b2[8-1] - 348.4 - 0.8) < 1e-2


    peptide.set_sequence('HYSHVDCPGHADYIK')
    peptide.create_fragmentation_pattern(R, bions=False, yions=False, xions=True)
    peptide.charge = 1
    ch = peptide.charge
    x =  [ ( pred + (ch -1)*R.mass_H)/ch for pred in peptide.x_series]

    assert abs(x[2] - 1380.6007 ) < 1e-4
    assert abs(x[3] - 1243.54179 ) < 1e-4
    assert abs(x[4] - 1144.47338) < 1e-4
    assert abs(x[-3] - 449.2400 ) < 1e-4
    assert abs(x[-4] - 564.267 ) < 1e-4

    peptide.set_sequence('GAGHSPHHVNGADTALQK')
    peptide.create_fragmentation_pattern(R, bions=True, yions=False, aions=False)
    print (peptide.b_series)
    return 1

    peptide.set_sequence('GAGHSPHHVNGADTALQK')
    peptide.create_fragmentation_pattern(R, bions=True, yions=False, aions=True)
    print (peptide.allseries)
    ch = 2
    charged =  [ ( pred + (ch -1)*R.mass_H)/ch for pred in peptide.allseries]
    print ("===")
    print (charged)

