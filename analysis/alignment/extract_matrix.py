#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
from sys import stdout, maxsize
import csv
import sys

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
        
def handle_args():
    usage="Generates matrix with flexible columns from featurealigner.tsv or featurealigner_requant.tsv file." \
          "For filtering or peak type highlighting use compute_full_matrix.py"
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('--in', dest="infile", required=True, help = 'feature aligner file')
    parser.add_argument("--out", required=True, default="", help="output matrix")
    parser.add_argument('--columns', default=['Intensity','RT','m_score'],  nargs = "+", help="Which columns are written")

    args = parser.parse_args(sys.argv[1:])

    return args

def main(options):
    import time

    # Read the files
    start = time.time()
    columns = options.columns
    print("Reading", options.infile,"extracting columns", columns)
    files = set()
    ids = set()
    entries = {}
    with open(options.infile) as f:
        f_csv = csv.reader(f,delimiter="\t")
        headers = next(f_csv)
        for row in f_csv:
            file = row[headers.index('align_origfilename')].split("/")[-1]
            files.add(file)
            id = row[headers.index('transition_group_id')] + "\t" + row[headers.index('ProteinName')]
            ids.add(id)
            entry = []
            for column in columns:
                entry.append(row[headers.index(column)])
            entries[id+file] = entry

    print("Writing", options.out)
    with open(options.out,'w') as f:
        #write header
        f.write("transition_group_id\tProteinName\t")
        for file in files:
            for column in columns:
                f.write("%s_%s\t"%(column,file))
        f.write("\n")

        #write content
        for id in ids:
            f.write(id+"\t")
            for file in files:
                if id+file in entries:
                    f.write("\t".join(entries[id+file])+"\t")
                else:
                    f.write("NA\t" * len(columns))
            f.write("\n")

    print("Took %ss" % (time.time() - start) )

if __name__=="__main__":
    options = handle_args()
    main(options)
