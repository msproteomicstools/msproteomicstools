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
$Maintainer: Lorenz Blum$
$Authors: Pedro Navarro, Lorenz Blum$
--------------------------------------------------------------------------
"""

from __future__ import print_function
import shutil
import sys
import getopt
import re


def pseudoreverse(seq, cleavage_rule):
    split = re.split("(%s)" % cleavage_rule, seq)
    return "".join([s[::-1] for s in split])


def usage():
    print("""Usage: python pseudoreverseDB [options]
Options:
-i   input.fasta input fasta file
-o   output.fasta output file with decoys
[-t] tag for decoys (optional, default 'DECOY_')
[-c] cleavage rule (optional, default '[KR](?!P)' =trypsin)
Creates pseudo-reversed decoys (peptide level reversion). If no cleavage site is found protein is simply reversed.
See http://dx.doi.org/10.1038/nmeth1019 Figure 6""")
    sys.exit(1)


def main(argv):
    # Get options
    try:
        opts, args = getopt.getopt(argv, "i:o:t:c:")
    except getopt.GetoptError:
        usage()

    decoytag = "DECOY_"
    cleavage_rule = '[KR](?!P)'
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
        if opt in ("-i", "--input"):
            i = arg
        if opt in ("-o", "-output"):
            o = arg
        if opt in ("-t", "--tag"):
            decoytag = arg
        if opt in ("-c", "--cleavage_rule"):
            cleavage_rule = arg

    # copy orig
    try:
        shutil.copy(i, o)
    except UnboundLocalError:
        usage()

    # append decoys
    out = open(o, 'a')
    seq = ""
    for line in open(i):
        if line.startswith('>'):
            out.write(pseudoreverse(seq, cleavage_rule) + '\n')
            seq = ""
            out.write('>%s%s' % (decoytag, line[1:]))
        else:
            seq += line.strip()
    #flush last one
    out.write(pseudoreverse(seq, cleavage_rule) + '\n')


if __name__ == '__main__':
    main(sys.argv[1:])



