#!/usr/bin/python
# -*- coding: utf-8  -*-

"""
 *
 * Program       : Select best peakgroup
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 13.03.2012
 *
 *
 * Copyright (C) 2011 Hannes Roest
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
 *
"""

import csv, sys
filename = sys.argv[1]

reader = csv.reader(open(filename), delimiter='\t')
header = reader.next()
header_d = dict( [ (h,i) for i,h in enumerate(header)] )

#print header
lines = list(reader)

# try to estimate an empirical FDR and only keep those lines above the fdr
decoys_ = 0
peptides = {}
proteins = {}
pep_per_prot = {}
for i,line in enumerate(lines):
    peptides[ line[header_d['Sequence'] ]] = 0
    proteins[ line[header_d['ProteinName'] ]] = 0
    if not pep_per_prot.has_key( line[header_d['ProteinName'] ]):
      pep_per_prot[  line[header_d['ProteinName'] ] ] = {}
    pep_per_prot[  line[header_d['ProteinName'] ] ][ line[header_d['Sequence'] ]] = 0

single_hits = len([0 for k,v in pep_per_prot.iteritems() if len(v) == 1])
print "Found precursors" , i, " and peptides", len(peptides), " and proteins", len(proteins), "(single hits", single_hits, "multiple hits", len(proteins) - single_hits, ")"




