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


def filterChromByTSV(infile, outfile, csvfile):

    import sqlite3
    conn = sqlite3.connect(infile)
    c = conn.cursor()

    import csv
    print("CSV file:", csvfile)
    print("outfile file:", outfile)
    rr = csv.reader(open(csvfile), delimiter="\t")
    header = rr.next()
    headerdict = dict([ (k,i) for i,k in enumerate(header) ])

    labels = []
    for line in rr:
        labels.extend( line[ headerdict['aggr_Fragment_Annotation']   ].split(";") )
        labels.extend( line[ headerdict['aggr_prec_Fragment_Annotation']   ].split(";") )

    labels = [ "'" + l + "'" for l in labels]
    labels_stmt = get_ids_stmt(labels)

    stmt = "SELECT ID FROM CHROMATOGRAM WHERE NATIVE_ID IN %s" % labels_stmt
    keep_ids = [i[0] for i in list(c.execute(stmt))]
    print("Keep %s chromatograms" % len(keep_ids) )

    nr_chrom = list(c.execute("SELECT COUNT(*) FROM CHROMATOGRAM"))[0][0]
    nr_spec = list(c.execute("SELECT COUNT(*) FROM SPECTRUM"))[0][0]

    assert(nr_chrom > 0)
    assert(nr_spec == 0)

    copyDatabase(c, conn, outfile, keep_ids)

def copyDatabase(c, conn, outfile, keep_ids):
    c.execute("ATTACH DATABASE '%s' AS other;" % outfile)

    # Tables: 
    #  - DATA
    #  - SPECTRUM
    #  - RUN
    #  - RUN_EXTRA
    #  - CHROMATOGRAM
    #  - PRODUCT
    #  - PRECURSOR
    copy_table(c, conn, keep_ids, "PRECURSOR", "CHROMATOGRAM_ID")
    copy_table(c, conn, keep_ids, "PRODUCT", "CHROMATOGRAM_ID")
    copy_table(c, conn, keep_ids, "DATA", "CHROMATOGRAM_ID")
    copy_table(c, conn, keep_ids, "CHROMATOGRAM", "ID")

    c.execute("CREATE TABLE other.RUN AS SELECT * FROM RUN");
    c.execute("CREATE TABLE other.SPECTRUM AS SELECT * FROM SPECTRUM");
    c.execute("CREATE TABLE other.RUN_EXTRA AS SELECT * FROM RUN_EXTRA");

    c.execute("CREATE INDEX other.data_chr_idx ON DATA(CHROMATOGRAM_ID);")
    c.execute("CREATE INDEX other.data_sp_idx ON DATA(SPECTRUM_ID);")
    c.execute("CREATE INDEX other.spec_rt_idx ON SPECTRUM(RETENTION_TIME);")
    c.execute("CREATE INDEX other.spec_mslevel ON SPECTRUM(MSLEVEL);")
    c.execute("CREATE INDEX other.spec_run ON SPECTRUM(RUN_ID);")
    c.execute("CREATE INDEX other.chrom_run ON CHROMATOGRAM(RUN_ID);")

    conn.commit()

def drop_table(c, conn, keep_ids, tbl, id_col):
    stmt = "CREATE TABLE %s_TMP AS SELECT * FROM %s WHERE %s IN " % (tbl, tbl, id_col)
    stmt += get_ids_stmt(keep_ids) + ";"
    c.execute(stmt)
    conn.commit()
    # 
    stmt = "DROP TABLE %s;" % tbl
    c.execute(stmt)
    conn.commit()
    # 
    stmt = "CREATE TABLE %s AS SELECT * FROM %s_TMP;" % (tbl, tbl)
    c.execute(stmt)
    conn.commit()
    #
    stmt = "DROP TABLE %s_TMP;" % tbl
    c.execute(stmt)
    conn.commit()

def copy_table(c, conn, keep_ids, tbl, id_col):
    stmt = "CREATE TABLE other.%s AS SELECT * FROM %s WHERE %s IN " % (tbl, tbl, id_col)
    stmt += get_ids_stmt(keep_ids) + ";"
    c.execute(stmt)
    conn.commit()

def get_ids_stmt(keep_ids):
    ids_stmt = "("
    for myid in keep_ids:
        ids_stmt += str(myid) + ","
    ids_stmt = ids_stmt[:-1]
    ids_stmt += ")"
    return ids_stmt 

def filterChromByNativeID(infile, outfile, selector):

    import sqlite3
    conn = sqlite3.connect(infile)
    c = conn.cursor()

    stmt = "SELECT ID FROM CHROMATOGRAM WHERE NATIVE_ID LIKE '%s%%'" % selector
    # keep_ids = list(c.execute(stmt))
    keep_ids = [i[0] for i in list(c.execute(stmt))]

    nr_chrom = list(c.execute("SELECT COUNT(*) FROM CHROMATOGRAM"))[0][0]
    nr_spec = list(c.execute("SELECT COUNT(*) FROM SPECTRUM"))[0][0]

    assert(nr_chrom > 0)
    assert(nr_spec == 0)

    copyDatabase(c, conn, outfile, keep_ids)

####### PARAMETERS
def handle_args():
    from optparse import OptionParser, OptionGroup
    import sys
    usage = 'A script to merge multiple sqMass files'

    import ast, argparse
    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", required=True, help='An input file', metavar="INP")
    parser.add_argument('--out', dest="outfile", required=True, help='An output file', metavar="OUT")
    parser.add_argument('--chrom_filter', dest="native_id_filter", help='Filter chromatograms by native id', metavar="NATIVE_ID")
    parser.add_argument('--tsv_filter', dest="tsv_filter", help='Filter chromatograms by TSV file', metavar="TSV")
    # parser.add_argument("--out", dest="outfile", required=True, default="", help="Output file")
    args = parser.parse_args(sys.argv[1:])

    return args

options = handle_args()

infile = options.infile

if options.native_id_filter:
    filterChromByNativeID(infile, options.outfile, options.native_id_filter)
elif options.tsv_filter:
    filterChromByTSV(infile, options.outfile, options.tsv_filter)

