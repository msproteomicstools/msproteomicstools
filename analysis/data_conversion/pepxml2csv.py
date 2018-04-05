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
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""
from __future__ import print_function

import sys
import os.path
#from msproteomicstoolslib.format import pepXMLReader
import csv

maxInt = sys.maxsize
decrement = True
while decrement:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.
    # http://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072

    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True


from pyteomics import pepxml

infile  = sys.argv[1]

outfile = os.path.splitext(infile)[0] + '.csv'

reader = pepxml.read(infile)

writer = csv.writer(open(outfile, 'w'), delimiter='\t')

## MYRIMATCH
{
    'end_scan': 1380,
    'retention_time_sec': 5190.16999999998,
    'index': 160,
    'assumed_charge': 2,
    'spectrum': '5P_HDMSE_121214_20.1380.1380.2',
    'search_hit': [
        {
            'hit_rank': 1,
            'calc_neutral_pep_mass': 3123.6467143833,
            'modifications': [],
            'modified_peptide': 'GKFDFSGNLDLEAFVLMAAEIGLWVILR',
            'peptide': 'GKFDFSGNLDLEAFVLMAAEIGLWVILR',
            'massdiff': 0.05699714826,
            'search_score': {
                'mvh': 12.1394600401,
                'number of matched peaks': 4.0,
                'xcorr': 0.9391081889259687,
                'number of unmatched peaks': 51.0,
                'mzFidelity': 10.015543238385
            },
            'proteins': [
                {
                    'num_tol_term': 2,
                    'protein': 'sp|Q8IW92|GLBL2_HUMAN',
                    'peptide_next_aa': 'P',
                    'protein_descr': None,
                    'peptide_prev_aa': 'R',
                },
                {
                    'protein': 'sp|Q8NCI6|GLBL3_HUMAN'
                }
            ],
            'num_missed_cleavages': 1,
            'tot_num_ions': 55,
            'num_tot_proteins': 2,
            'num_matched_ions': 4
        }
    ],
    'precursor_neutral_mass': 3123.58971723504,
    'num_decoy_comparisons': '0',
    'start_scan': 1380,
    'num_target_comparisons': '936',
    'spectrumNativeID': 'controllerType=0 controllerNumber=1 scan=1380'
}


## X!Tandem

{'end_scan': 113,
'retention_time_sec': 963.599,
'index': 367,
'assumed_charge': 2,
'spectrum': '5P_HDMSE_121214_20.00113.00113.2',
'search_hit':
    [{
    'hit_rank': 1,
    'calc_neutral_pep_mass': 3139.3787,
    'modifications': [{'position': 20, 'mass': 160.0306}],
    'modified_peptide': 'STDVDAVPYTGADSTQGTWCEDEPVGARR',
    'peptide': 'STDVDAVPYTGADSTQGTWCEDEPVGARR',
    'num_matched_ions': 1,
    'search_score': {
                    'zscore': 0.0,
                    'xscore': 0.0,
                    'bscore': 1.0,
                    'yscore': 0.0,
                    'cscore': 0.0,
                    'hyperscore': 229.0,
                    'ascore': 0.0,
                    'expect': 53.0,
                    'nextscore': 277.0},
    'proteins': [{
                'num_tol_term': 2,
                'protein': 'DECOY_Q49AJ0',
                'peptide_next_aa': 'V',
                'protein_descr': None,
                'peptide_prev_aa': 'K'}],
    'num_missed_cleavages': 0,
    'tot_num_ions': 56,
    'num_tot_proteins': 1,
    'is_rejected': False,
    'massdiff': -0.03}],
'precursor_neutral_mass': 3139.349,
'start_scan': 113}

## OMSSA
{
'end_scan': 1380,
'index': 92,
'assumed_charge': 2,
'spectrum': '5P_HDMSE_121214_20.01380.01380.2',
'search_hit': [{
            'hit_rank': 1,
            'calc_neutral_pep_mass': 3125.59204101562,
            'modifications': [],
            'modified_peptide':
            'DQALSISAMEELPRVLYLPLFMEAFSR',
            'peptide': 'DQALSISAMEELPRVLYLPLFMEAFSR',
            'num_matched_ions': 2,
            'search_score': {
                            'pvalue': 129.44166834226553,
                            'expect': 3883.250050267966},
            'proteins': [{
                        'num_tol_term': 0,
                        'protein': 'O95521',
                        'peptide_next_aa': 'R',
                        'protein_descr': 'PRAME family member 1 OS=Homo sapiens GN=PRAMEF1 PE=2 SV=2',
                        'peptide_prev_aa': 'R'}],
            'tot_num_ions': 52,
            'num_tot_proteins': 1,
            'is_rejected': False,
            'massdiff': 0.014892578125}],
'precursor_neutral_mass': 3125.60693359375,
'start_scan': 1380}


{
'end_scan': 1380,
'index': 92,
'assumed_charge': 2,
'spectrum': '5P_HDMSE_121214_20.01380.01380.2',
'search_hit': [{
            'hit_rank': 1,
            'calc_neutral_pep_mass': 3125.59204101562,
            'modifications': [],
            'modified_peptide':
            'DQALSISAMEELPRVLYLPLFMEAFSR',
            'peptide': 'DQALSISAMEELPRVLYLPLFMEAFSR',
            'num_matched_ions': 2,
            'search_score': {
                            'pvalue': 129.44166834226553,
                            'expect': 3883.250050267966},
            'proteins': [{
                        'num_tol_term': 0,
                        'protein': 'O95521',
                        'peptide_next_aa': 'R',
                        'protein_descr': 'PRAME family member 1 OS=Homo sapiens GN=PRAMEF1 PE=2 SV=2',
                        'peptide_prev_aa': 'R'}],
            'tot_num_ions': 52,
            'num_tot_proteins': 1,
            'is_rejected': False,
            'massdiff': 0.014892578125}],
'precursor_neutral_mass': 3125.60693359375,
'start_scan': 1380
}



modifications_example = [{'position': 20, 'mass': 160.0306}]
iterable_items = ['search_score']
iterable_listdict_items = ['proteins']
header = ['start_scan','end_scan','index','assumed_charge','spectrum','precursor_neutral_mass']

headerset = False
headerdict = False

for hit in reader :

    if not 'search_hit' in hit : continue

    if not headerset : #read the elements the header is going to contain
        non_iterable_search_hit = [ i for i in hit['search_hit'][0].keys() ]
        for el in iterable_items :
            if el in non_iterable_search_hit : non_iterable_search_hit.remove(el)
            non_iterable_search_hit.extend(hit['search_hit'][0][el])
        for el in iterable_listdict_items :
            if el in non_iterable_search_hit : non_iterable_search_hit.remove(el)
            non_iterable_search_hit.extend([i for i in hit['search_hit'][0][el][0].keys()])
        header.extend(non_iterable_search_hit)
        #Generate a header dictionary
        header_dict = dict([ (l,i) for i,l in enumerate(header)])
        writer.writerow(header)
        headerset = True

    new_row = []
    for element in hit :
        if element in header_dict : new_row.append( (header_dict[element] , hit[element]) )
    for id, element in hit['search_hit'][0].items() :
        if id == 'modifications' :
            #Pass the list of dictionaries as a string
            mods = ','.join(str(x) for x in element)
            if id in header_dict : new_row.append( (header_dict[id], mods) )
            continue
        if id in header_dict : new_row.append( (header_dict[id], element) )
    for itelement in iterable_items :
        for id,element in hit['search_hit'][0][itelement].items() :
            if id in header_dict : new_row.append( (header_dict[id], element) )
    for itelement in iterable_listdict_items :
        dict_elements = {}
        for d in hit['search_hit'][0][itelement] :
            for id,el in d.items() :
                if id in dict_elements : dict_elements[id] += "#" + str(el)
                else : dict_elements[id] = str(el)
        for id, element in dict_elements.items() :
            if id in header_dict : new_row.append( (header_dict[id], element) )


    new_row.sort(key=lambda x: x[0])
    new_row_formatted = [ y for (x,y) in new_row ]

    #Check the final size of the row
    if len(new_row_formatted) < len(header) :
        new_row_formatted.extend(['-' for x in range(len(header)-len(new_row_formatted))])

    writer.writerow(new_row_formatted)

print("...file %s written." % outfile)
