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
$Maintainer: Hannes Roest $
$Authors: Hannes Roest $
--------------------------------------------------------------------------
"""

"""
Script to filter peakview csv libraries

Applies to PeakView libraries that are unfiltered and have not gone through the
conversion process with the spectrast2tsv.py

- remove transitions that are not within the cutoff (less than 350 m/z)
- normalize intensity to 10k 
- fix modifications

    python filterPeakview.py input.csv output.csv


"""

import sys,csv

inp = sys.argv[1]
outp = sys.argv[2]

nr_removed_seqs = 0

r = csv.reader(open(inp), delimiter="\t")
wr = csv.writer(open(outp, "w"), delimiter="\t")

header = next(r)
# print(header)
hdict = dict( [ (c,i) for i,c in enumerate(header) ] )
wr.writerow( header + ["UniProtID"] )

def filterNormalizeTransitions(stack, topN):
    """
    Process a single group of transitions and filter them as follows:

        - remove transitions within +/- 10 Da of the precursor
        - remove transitions below 350 m/z

    It then selects the topN transitions from the remaining transitions.

    It also normalizes the intensity to 10 000 and the uniprot id from
    "sp|P04114|APOB_HUMAN" to "P04114".

    """
    try:
        stack.sort(lambda x,y: -cmp( float(x[hdict["relative_intensity"]]),  float(y[hdict["relative_intensity"]]) ) )
    except ValueError:
        print "Value Error in sort, abort for", stack[0]
        return []

    # Normalize uniprot id 
    prec = float(stack[0][hdict["Q1"]] )
    unipr = stack[0][hdict["uniprot_id"]]
    if unipr.startswith("sp"):
        unipr = unipr.split("sp|")[1].split("|")[0]

    # Filter transitions by m/z and around precursor
    MZCUTOFF = 350.0
    stack_ = [s for s in stack if float(s[hdict["Q3"]]) > MZCUTOFF and 
            (float(s[hdict["Q3"]])  < prec - 10 or float(s[hdict["Q3"]]) > prec +10)]

    if len(stack_) == 0:
        print "abort for" , stack[0][hdict["modification_sequence"]]
        # print stack
        return []

    stack = stack_

    if len(stack) == 0:
        return []

    maxint = float( stack[0][hdict["relative_intensity"]] )
    # select top N
    stack = stack[:topN]

    # normalize intensity and append uniprot id
    for ll in stack: 
        ll[hdict["relative_intensity"]] = 10000 * float(ll[hdict["relative_intensity"]]) / maxint
        ll.append(unipr)

    return stack

missingMods = set([])

def mapModifications(stack):
    """
    Map modifications from Peakview format to OpenSWATH format.

    Currently not supported modifications:

    '[Ntr]', '[Myr]', '[Pty]', '[PGQ]', '[3Ox]', '[hHm]', '[AAS]', '[Orn]',
    '[GCn]', '[1Br]', '[Gsl]', '[Dhy]', '[Pup]', '[HKy]', '[MeP]', '[Trx]',
    '[K1G]', '[2HF]', '[Kyn]', '[NGc]', '[Oxi]', '[SU1]', '[GPE]', '[Oct]',
    '[G2Z]', '[PGP]', '[1G2]', '[Hex]', '[MSH]', '[Pho]', '[Pgl]', '[SuO]',
    '[CEt]', '[MDe]', '[Iod]', '[Poy]', '[AGA]', '[Hep]', '[DTM]', '[NHH]',
    '[Frm]', '[2Br]', '[Btn]', '[CHD]', '[HNE]', '[G0c]', '[6HF]', '[G2Y]',
    '[-1K]', '[KXX]', '[Sud]', '[FMN]', '[Hps]', '[PyP]', '[5Hx]', '[GRL]',
    '[AAA]', '[GMP]', '[cGP]', '[K2P]', '[SU2]', '[CRM]', '[2NF]', '[GPh]',
    '[Lip]', '[1K1]', '[NHN]', '[HHN]', '[K2F]', '[h3M]', '[AAR]', '[NaX]',
    '[DPO]', '[UMP]', '[Adh]', '[5HF]', '[h1M]', '[K1F]', '[2Cl]', '[NAc]',
    '[ddR]', '[G0X]', '[3Me]', '[pHc]', '[dOx]', '[2CM]', '[GGQ]', '[G1F]',
    '[Ret]', '[G2W]', '[Byr]', '[-1R]', '[G0N]', '[2Me]', '[1G1]', '[K3e]',
    '[Dec]', '[Qin]', '[AMP]', '[Ahp]', '[NF2]', '[Suc]', '[dAm]', '[4HF]',
    '[PGE]', '[G2d]', '[G1c]', '[Pal]', '[K1H]', '[Hmn]', '[1K3]', '[OH2]',
    '[SSe]', '[1Ac]', '[K1P]', '[Ami]', '[YDA]', '[G0F]', '[UGG]', '[CuX]',
    '[dHx]', '[Cro]', '[pHx]', '[PAI]', '[3HF]', '[1Me]', '[Cox]', '[NFX]',
    '[G1S]', '[2Hx]', '[OH1]', '[Dea]', '[-2H]', '[1Cl]', '[+1K]', '[HNc]',
    '[1HF]', '[PYD]', '[Amn]', '[6Hx]', '[1HN]', '[GdP]', '[Hse]', '[2Ox]',
    '[PRC]', '[+1R]', '[pRb]', '[MPr]', '[Dip]', '[2M+]', '[K1e]', '[LAA]',
    '[3Hx]', '[hCo]', '[CAM]', '[HEM]'

    """

    global nr_removed_seqs 
    global missingMods 
    if len(stack) == 0:
        return stack

    # mcopy = stack[0][hdict["modification_sequence"]][:]
    mseq = stack[0][hdict["modification_sequence"]]
    if mseq.find("[") == -1:
        return stack

    mymods = [
        ['C[160]',   4,  '[CAM]'],
        ['M[147]',  35,  '[Oxi]'],
        ['W[202]',  35,  '[Oxi]'],
        ['H[153]',  35,  '[Oxi]'],
        ['K[136]', 259,  '[+08]'],
        ['R[166]', 267,  '[+10]'],
        ['E[111]',  27,  '[PGE]'],
        ['Q[111]',  28,  '[PGQ]'],
        ['C[143]',  26,  '[PCm]'],
        ['n[43]',    5,  '[CRM]'],
        ['S[167]',  21,  '[Pho]'],
        ['T[181]',  21,  '[Pho]'],
        ['Y[243]',  21,  '[Pho]'],
        ['N[115]',   7,  '[Dea]'],
        ['Q[129]',   7,  '[Dea]'],
        ['C[149]',  39,  '[XXX]'],
        ['D[131]',  35,  '[Oxi]'],
        ['K[144]',  35,  '[Oxi]'],
        ['Y[179]',  35,  '[Oxi]'],
        ['R[172]',  35,  '[Oxi]'],
        ['N[130]',  35,  '[Oxi]'],
        ['P[113]',  35,  '[Oxi]'],
        ['C[119]',  35,  '[Oxi]'] ]

    for a,b,c in mymods:
        aa = a[0]
        mseq = mseq.replace( aa+c, aa + "(UniMod:%s)" % b)

    for ll in stack: 
        ll[hdict["modification_sequence"]] = mseq


    import re
    allHits = re.compile("\[[^\]]*\]")
    expr = re.compile("\[[^\]]*([^0-9\]]+)[^\]]*\]")
    if expr.search(mseq):
        nr_removed_seqs += 1
        ## print "before", mcopy, "after", mseq, "remove seq"
        ## print allHits.findall(mseq)
        missingMods.update( allHits.findall(mseq) )
        return []

    return stack


###################################
## parse through all data and process one transition group at a time
stack = []
prev_unique_id = ""
for l in r:
    unique_id = l[hdict["modification_sequence"]] + "/" + l[hdict["prec_z"]]

    if prev_unique_id == "":
        prev_unique_id = unique_id
    if unique_id != prev_unique_id:
        st = filterNormalizeTransitions(stack, 6)
        st = mapModifications(st)
        wr.writerows(st)
        prev_unique_id = unique_id
        stack = []

    stack.append(l)
    

# Last set
st = filterNormalizeTransitions(stack, 6)
wr.writerows(st)

print "removed sequences", nr_removed_seqs 
# print missingMods

