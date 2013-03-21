#!/usr/bin/env python
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
filename_out = sys.argv[2]
sort_by = sys.argv[3]
decoy_name = sys.argv[4]
fdr_cutoff = float(sys.argv[5])
reverse = sys.argv[6]
reverse = (reverse == 'TRUE')
use_fdr_score = sys.argv[7]
use_fdr_score = (use_fdr_score == 'TRUE')

reader = csv.reader(open(filename), delimiter='\t')
header = reader.next()
header_d = dict( [ (h,i) for i,h in enumerate(header)] )
try:
  sort_pos = header_d[sort_by]
except KeyError:
    print "Could not find key", sort_by, "in", header_d
    sys.exit()


print "Sort by '%s' (position %s) and reverse %s" % (sort_by, sort_pos, reverse)
lines = list(reader)
if reverse:
  lines.sort( lambda x,y: -cmp( float(x[sort_pos]), float(y[sort_pos]) ) )
else:
  lines.sort( lambda x,y: cmp( float(x[sort_pos]), float(y[sort_pos]) ) )

# try to estimate an empirical FDR and only keep those lines above the fdr
decoys_ = 0
for i,line in enumerate(lines):
    is_decoy = False
    for cell in line: 
        if cell.find(decoy_name) != -1: 
            is_decoy = True

    if is_decoy:
        decoys_ += 1

    if use_fdr_score:
      fdr = float(line[sort_pos])
    elif (i+1-decoys_) > 0:
      fdr = decoys_ * 1.0/ (i+1-decoys_)
    else: fdr = 1.0
    line.append(fdr)
    if(fdr > fdr_cutoff): break


print "Stop after %s entries of which %s are decoys (est. fdr=%s)" % (i, decoys_, fdr)
above_cutoff = lines[:i]

header.append('FDREstimate')
writer = csv.writer(open(filename_out, 'w'), delimiter='\t')
writer.writerow(header)
writer.writerows(above_cutoff)


