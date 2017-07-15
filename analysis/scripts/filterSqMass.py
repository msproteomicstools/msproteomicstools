

def filterChromByTSV(infile, csvfile):

    import sqlite3
    conn = sqlite3.connect(infile)
    c = conn.cursor()

    import csv
    print("CSV file:", csvfile)
    rr = csv.reader(open(csvfile), delimiter="\t")
    header = rr.next()
    headerdict = dict([ (k,i) for i,k in enumerate(header) ])

    labels = []
    for line in rr:
        labels.extend( line[ headerdict['aggr_Fragment_Annotation']   ].split(";") )

    labels = [ "'" + l + "'" for l in labels]
    labels_stmt = get_ids_stmt(labels)

    stmt = "SELECT ID FROM CHROMATOGRAM WHERE NATIVE_ID IN %s" % labels_stmt
    keep_ids = [i[0] for i in list(c.execute(stmt))]
    print("Keep %s chromatograms" % len(keep_ids) )

    nr_chrom = list(c.execute("SELECT COUNT(*) FROM CHROMATOGRAM"))[0][0]
    nr_spec = list(c.execute("SELECT COUNT(*) FROM SPECTRUM"))[0][0]

    assert(nr_chrom > 0)
    assert(nr_spec == 0)

    drop_table(c, conn, keep_ids, "PRECURSOR", "CHROMATOGRAM_ID")
    drop_table(c, conn, keep_ids, "PRODUCT", "CHROMATOGRAM_ID")
    drop_table(c, conn, keep_ids, "DATA", "CHROMATOGRAM_ID")
    drop_table(c, conn, keep_ids, "CHROMATOGRAM", "ID")

    c.execute("VACUUM;")
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


def get_ids_stmt(keep_ids):
    ids_stmt = "("
    for myid in keep_ids:
        ids_stmt += str(myid) + ","
    ids_stmt = ids_stmt[:-1]
    ids_stmt += ")"
    return ids_stmt 

def filterChromByNativeID(infile, selector):

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

    drop_table(c, conn, keep_ids, "PRECURSOR", "CHROMATOGRAM_ID")
    drop_table(c, conn, keep_ids, "PRODUCT", "CHROMATOGRAM_ID")
    drop_table(c, conn, keep_ids, "DATA", "CHROMATOGRAM_ID")
    drop_table(c, conn, keep_ids, "CHROMATOGRAM", "ID")

    c.execute("VACUUM;")
    conn.commit()


####### PARAMETERS
def handle_args():
    from optparse import OptionParser, OptionGroup
    import sys
    usage = 'A script to merge multiple sqMass files'

    import ast, argparse
    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", required=True, help='An input file', metavar="INP")
    parser.add_argument('--chrom_filter', dest="native_id_filter", help='Filter chromatograms by native id', metavar="NATIVE_ID")
    parser.add_argument('--tsv_filter', dest="tsv_filter", help='Filter chromatograms by TSV file', metavar="TSV")
    # parser.add_argument("--out", dest="outfile", required=True, default="", help="Output file")
    args = parser.parse_args(sys.argv[1:])

    return args

options = handle_args()

infile = options.infile

print options.native_id_filter

if options.native_id_filter:
    filterChromByNativeID(infile, options.native_id_filter)
elif options.tsv_filter:
    filterChromByTSV(infile, options.tsv_filter)

