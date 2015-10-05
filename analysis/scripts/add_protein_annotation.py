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

import sys
import csv

"""
Doc :
    A script to annotate a csv with protein descriptions from FASTA.
    It will search for a column named ProteinName, look it up
    in the FASTA file and add a column with the name ProteinDescription to the csv.
"""

####### PARAMETERS
def handle_args():
    from optparse import OptionParser, OptionGroup
    usage = 'A script to annotate a csv with protein descriptions from FASTA. ' + \
            'It will search for a column named ProteinName, look it up' + \
            'in the FASTA file and add a column with the name ProteinDescription to the csv.'

    parser = OptionParser(usage=usage)
    group = OptionGroup(parser, "Feature Alignments Options", "")
    group.add_option("--fasta", dest="fasta_file", help="A fasta file", type="string")
    group.add_option("--in", dest="mprophet_file", help="A mProphet output file", type="string")
    group.add_option("--out", dest="outfile", help="A modified mProphet output file", type="string")

    parser.add_option_group(group)
    options, args = parser.parse_args(sys.argv[1:])
    return options

options = handle_args()

mpr = open(options.mprophet_file)
out = open(options.outfile, "w")

def map_fasta_identifier(fasta_file):
    # read in the mapping of id to description from FASTA
    mapping = {}
    from Bio import SeqIO
    handle = open(fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        mapping[record.id] = record.description

    handle.close()
    return mapping

mapping = map_fasta_identifier(options.fasta_file)
r = csv.reader(mpr, delimiter="\t")
we = csv.writer(out, delimiter="\t")
header = next(r)
header.append("ProteinDescription")
header_d = dict([(c,i) for i,c in enumerate(header)])
we.writerow(header)

for line in r:
    try:
        descr = mapping[ line[ header_d["ProteinName"] ] ] 
        line.append(descr)
    except KeyError:
        # sometimes there are RT peptides in there
        line.append("NA")
    we.writerow(line)


out.close()
mpr.close()

