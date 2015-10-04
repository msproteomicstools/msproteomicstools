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

####### PARAMETERS
def handle_args():
    from optparse import OptionParser, OptionGroup
    usage = """
    A script to display which proteins were not present in a library.
    It will search for a column named ProteinName in a csv, look it up
    in the FASTA file and then print out all extra entries in the FASTA file
    (except those starting with DECOY and sp|).
    """
    parser = OptionParser(usage=usage)
    group = OptionGroup(parser, "Feature Alignments Options", "")
    group.add_option("--fasta", dest="fasta_file", help="A fasta file", type="string")
    group.add_option("--in", dest="csv_file", help="A csv file with the column ProteinName", type="string")
    group.add_option("--out", dest="outfile", help="A csv output file containing all proteins not in the csv file", type="string")

    parser.add_option_group(group)
    options, args = parser.parse_args(sys.argv[1:])
    return options

options = handle_args()

def map_fasta_identifier(fasta_file):
    """
    Read in the mapping of id to description from FASTA
    """
    mapping = {}
    from Bio import SeqIO
    handle = open(fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        mapping[record.id] = record.description

    handle.close()
    return mapping

# Open the csv file, get the mapping
mpr = open(options.csv_file)
mapping = map_fasta_identifier(options.fasta_file)
r = csv.reader(mpr, delimiter="\t")
header = next(r)
header_d = dict([(c,i) for i,c in enumerate(header)])

# Collect all the proteins that were found
found = set([])
for line in r:
    found.update( [line[ header_d["ProteinName"] ] ])

# Compute the proteins that were not found
notfound = set(mapping.keys()).difference(found)
notfound = [n for n in notfound if not n.startswith("DECOY") and not n.startswith("sp|")]

# Write out the result
out = open(options.outfile, "w")
we = csv.writer(out, delimiter="\t")
for line in notfound:
    we.writerow( [line, mapping[line] ] )

out.close()


