
# convert a TSV output from feature alignment tsv (or mprophet output) 
# to Mayu input.

# TODO add modifications!

import sys, csv, argparse

import gzip


def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program converts an mProphet TSV output file to a Mayu input file"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", required=True, nargs = '+', help = 'A list of mProphet output files containing all peakgroups (use quotes around the filenames)')
    parser.add_argument("--out", dest="outfile", required=True, help="Output file")

    args = parser.parse_args(sys.argv[1:])
    return args

def get_filehandler(inputf):
    if inputf.endswith('.gz'):
        import gzip
        return gzip.open(inputf,'rb')
    else:
        return open(inputf)

args = handle_args()

outfile = csv.writer(open(args.outfile, "w"), delimiter=",")

if len(args.infiles) != 1:
    raise Exception("More than one infile not yet implemented!") 

# Input
inputf = args.infiles[0]
filehandler = get_filehandler(inputf)
infile = csv.reader(filehandler, delimiter="\t")

# Inputheader
header = infile.next()
header_dict = dict([(h, i) for i,h in enumerate(header)])

outfile.writerow(["Identifier", "Sequence", "Protein", "Mod", "MScore"])
line_nr = 0
for line in infile:
    line_nr += 1
    runid = line[ header_dict["align_runid"] ]
    #trgroup = line[ header_dict["transition_group_record"]  ]
    firstline = "run%s.%s.%s.%s" % (runid, line_nr, line_nr, line[ header_dict["Charge"]  ])
    newline = [ 
        firstline,
        line[ header_dict["Sequence"]  ],
        line[ header_dict["ProteinName"]  ],
        "", # no mods
        1 - float(line[ header_dict["m_score"]  ]),
    ]
    outfile.writerow(newline)

