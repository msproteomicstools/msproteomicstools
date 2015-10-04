
# convert a TSV output from feature alignment tsv to Mayu input.

# TODO add modifications!

import sys, csv, argparse

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program converts an mProphet TSV output file to a Mayu input file"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", required=True, help = 'A list of mProphet output files containing all peakgroups (use quotes around the filenames)')
    parser.add_argument("--out", dest="outfile", required=True, help="Output file")
    # parser.add_argument('--use_alignment_id', action='store_true', default=False, help="Use the column 'align_runid' to dermine from which run an identification originated")

    args = parser.parse_args(sys.argv[1:])
    return args

def get_filehandler(inputf):
    if inputf.endswith('.gz'):
        import gzip
        return gzip.open(inputf,'rb')
    else:
        return open(inputf)

def write_output(infile, outfile):
    for line_nr, line in enumerate(infile):
        runid = line[ header_dict["align_runid"] ]
        firstline = "run%s.%s.%s.%s" % (runid, line_nr, line_nr, line[ header_dict["Charge"]  ])
        newline = [ 
            firstline,
            line[ header_dict["Sequence"]  ],
            line[ header_dict["ProteinName"]  ],
            "", # no mods
            1 - float(line[ header_dict["m_score"]  ]),
        ]
        outfile.writerow(newline)

if __name__=="__main__":
    args = handle_args()

    outfile = csv.writer(open(args.outfile, "w"), delimiter=",")
    outfile.writerow(["Identifier", "Sequence", "Protein", "Mod", "MScore"])

    # Input
    inputf = args.infile
    filehandler = get_filehandler(inputf)
    infile = csv.reader(filehandler, delimiter="\t")

    # Inputheader
    header = next(infile)
    header_dict = dict([(h, i) for i,h in enumerate(header)])

    write_output(infile, outfile)

