

def mergeDB(basedb, otherdb):
    import sqlite3
    conn = sqlite3.connect(basedb)
    c = conn.cursor()

    newrun_start_ = list(c.execute("SELECT MAX(ID) FROM RUN"))[0][0]
    newspec_start_ = list(c.execute("SELECT MAX(ID) FROM SPECTRUM"))[0][0]
    newchrom_start_ = list(c.execute("SELECT MAX(ID) FROM CHROMATOGRAM"))[0][0]

    newrun_start = 0
    newspec_start = 0
    newchrom_start = 0
    if newrun_start_ is not None: 
        newrun_start = newrun_start_ + 1
    if newspec_start_ is not None: 
        newspec_start = newspec_start_ + 1
    if newchrom_start_ is not None: 
        newchrom_start = newchrom_start_ + 1


    print "attach", otherdb
    print "newspec start", newspec_start
    print "newchrom start", newchrom_start
    # NOTE: this currently assumes that the 
    # c.execute("BEGIN TRANSACTION")
    c.execute("ATTACH DATABASE '%s' AS other;" % otherdb)
    other_chrom = list(c.execute("SELECT COUNT(*) FROM CHROMATOGRAM"))[0][0]
    other_spec = list(c.execute("SELECT COUNT(*) FROM SPECTRUM"))[0][0]
    ###################################
    # do RUN
    stmt = "INSERT INTO RUN (ID, FILENAME, NATIVE_ID) SELECT ID + %s, FILENAME, NATIVE_ID FROM other.RUN;" % (newrun_start)
    c.execute(stmt)
    ###################################
    # do RUN_EXTRA
    stmt = "INSERT INTO RUN_EXTRA (RUN_ID, DATA) SELECT RUN_ID + %s, DATA FROM other.RUN_EXTRA;" % (newrun_start)
    c.execute(stmt)
    ###################################
    # do DATA
    stmt = """
    INSERT INTO DATA (SPECTRUM_ID, COMPRESSION, DATA_TYPE, DATA)
    SELECT SPECTRUM_ID + %s, COMPRESSION, DATA_TYPE, DATA 
    FROM other.DATA WHERE SPECTRUM_ID NOT NULL;
    """ % (newspec_start)
    if other_spec: c.execute(stmt)
    stmt = """
    INSERT INTO DATA (CHROMATOGRAM_ID, COMPRESSION, DATA_TYPE, DATA)
    SELECT CHROMATOGRAM_ID + %s, COMPRESSION, DATA_TYPE, DATA
    FROM other.DATA WHERE CHROMATOGRAM_ID NOT NULL;
    """ % (newchrom_start)
    if other_chrom: c.execute(stmt)
    ###################################
    # do PRECURSOR
    stmt = """
    INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, PEPTIDE_SEQUENCE, DRIFT_TIME, ACTIVATION_METHOD, ACTIVATION_ENERGY, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER)
    SELECT SPECTRUM_ID + %s, CHARGE, PEPTIDE_SEQUENCE, DRIFT_TIME, ACTIVATION_METHOD, ACTIVATION_ENERGY, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER
    FROM other.PRECURSOR WHERE SPECTRUM_ID NOT NULL;
    """ % (newspec_start)
    if other_spec: c.execute(stmt)
    stmt = """
    INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, PEPTIDE_SEQUENCE, DRIFT_TIME, ACTIVATION_METHOD, ACTIVATION_ENERGY, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER)
    SELECT CHROMATOGRAM_ID + %s, CHARGE, PEPTIDE_SEQUENCE, DRIFT_TIME, ACTIVATION_METHOD, ACTIVATION_ENERGY, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER
    FROM other.PRECURSOR WHERE CHROMATOGRAM_ID NOT NULL;
    """ % (newchrom_start)
    if other_chrom: c.execute(stmt)
    ###################################
    # do PRODUCT
    stmt = """
    INSERT INTO PRODUCT (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER)
    SELECT SPECTRUM_ID + %s, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER
    FROM other.PRODUCT WHERE SPECTRUM_ID NOT NULL;
    """ % (newspec_start)
    if other_spec: c.execute(stmt)
    stmt = """
    INSERT INTO PRODUCT (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER)
    SELECT CHROMATOGRAM_ID + %s, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER 
    FROM other.PRODUCT WHERE CHROMATOGRAM_ID NOT NULL;
    """ % (newchrom_start)
    if other_chrom: c.execute(stmt)
    ###################################
    # do SPECTRUM
    stmt = """
    INSERT INTO SPECTRUM (ID, RUN_ID, MSLEVEL, RETENTION_TIME, SCAN_POLARITY, NATIVE_ID) 
    SELECT ID + %s, RUN_ID + %s, MSLEVEL, RETENTION_TIME, SCAN_POLARITY, NATIVE_ID 
    FROM other.SPECTRUM;
    """ % (newspec_start, newrun_start)
    c.execute(stmt)
    ###################################
    # do CHROMATOGRAM
    stmt = """
    INSERT INTO CHROMATOGRAM (ID, RUN_ID, NATIVE_ID) 
    SELECT ID + %s, RUN_ID + %s, NATIVE_ID 
    FROM other.CHROMATOGRAM;
    """ % (newchrom_start, newrun_start)
    if other_chrom: c.execute(stmt)

    ## done
    conn.commit()
    #c.execute("END TRANSACTION")



####### PARAMETERS
def handle_args():
    from optparse import OptionParser, OptionGroup
    import sys
    usage = 'A script to merge multiple sqMass files'

    import ast, argparse
    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", required=True, nargs = '+', help = 'A list of files to be merged', metavar="INP")
    parser.add_argument("--out", dest="outfile", required=True, default="", help="Output file")
    args = parser.parse_args(sys.argv[1:])

    return args

options = handle_args()

infiles = options.infiles
if len(infiles) < 2:
    raise Exception("Need at least 2 input files")

from shutil import copyfile
copyfile(infiles[0], options.outfile)

infiles = infiles[1:]

for f in infiles:
    mergeDB(options.outfile, f)

# INSERT INTO RUN (ID, FILENAME, NATIVE_ID)  SELECT ID + (SELECT MAX(ID) FROM RUN) +1, FILENAME, NATIVE_ID FROM other.RUN;

"""
    void MzMLSqliteHandler::createTables()
    {
      // delete file if present
      QFile file (filename_.toQString());
      file.remove();

      sqlite3 *db = openDB();

      // Create SQL structure
      char const *create_sql =

        // data table
        //  - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        //  - data_type is one of 0 = mz, 1 = int, 2 = rt
        //  - data contains the raw (blob) data for a single data array
        "CREATE TABLE DATA(" \
        "SPECTRUM_ID INT," \
        "CHROMATOGRAM_ID INT," \
        "COMPRESSION INT," \
        "DATA_TYPE INT," \
        "DATA BLOB NOT NULL" \
        ");" \
        // CREATE INDEX data_chr_idx ON DATA(CHROMATOGRAM_ID);
        // CREATE INDEX data_sp_idx ON DATA(SPECTRUM_ID);

        // spectrum table
        "CREATE TABLE SPECTRUM(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "RUN_ID INT," \
        "MSLEVEL INT NULL," \
        "RETENTION_TIME REAL NULL," \
        "SCAN_POLARITY INT NULL," \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \
        // CREATE INDEX spec_rt_idx ON SPECTRUM(RETENTION_TIME);
        // CREATE INDEX spec_mslevel ON SPECTRUM(MSLEVEL);
        // CREATE INDEX spec_run ON SPECTRUM(RUN_ID);

        // ms-run table
        "CREATE TABLE RUN(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "FILENAME TEXT NOT NULL, " \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \

        // ms-run extra table
        "CREATE TABLE RUN_EXTRA(" \
        "RUN_ID INT," \
        "DATA BLOB NOT NULL" \
        ");" \
        // CREATE INDEX run_extra ON RUN_EXTRA(RUN_ID);

        // chromatogram table
        "CREATE TABLE CHROMATOGRAM(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "RUN_ID INT," \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \
        // CREATE INDEX chrom_run ON CHROMATOGRAM(RUN_ID);



        // product table
        "CREATE TABLE PRODUCT(" \
        "SPECTRUM_ID INT," \
        "CHROMATOGRAM_ID INT," \
        "CHARGE INT NULL," \
        "ISOLATION_TARGET REAL NULL," \
        "ISOLATION_LOWER REAL NULL," \
        "ISOLATION_UPPER REAL NULL" \
        ");" \
        // CREATE INDEX product_chr_idx ON DATA(CHROMATOGRAM_ID);
        // CREATE INDEX product_sp_idx ON DATA(SPECTRUM_ID);

        // precursor table
        "CREATE TABLE PRECURSOR(" \
        "SPECTRUM_ID INT," \
        "CHROMATOGRAM_ID INT," \
        "CHARGE INT NULL," \
        "PEPTIDE_SEQUENCE TEXT NULL," \
        "DRIFT_TIME REAL NULL," \
        "ACTIVATION_METHOD INT NULL," \
        "ACTIVATION_ENERGY REAL NULL," \
        "ISOLATION_TARGET REAL NULL," \
        "ISOLATION_LOWER REAL NULL," \
        "ISOLATION_UPPER REAL NULL" \
        ");";
        // CREATE INDEX precursor_chr_idx ON DATA(CHROMATOGRAM_ID);
        // CREATE INDEX precursor_sp_idx ON DATA(SPECTRUM_ID);

      // Execute SQL statement
      char *zErrMsg = 0;
      int rc;
      rc = sqlite3_exec(db, create_sql, callback, 0, &zErrMsg);
      if (rc != SQLITE_OK)
      {
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
            zErrMsg);
      }
      sqlite3_close(db);
    }
"""
