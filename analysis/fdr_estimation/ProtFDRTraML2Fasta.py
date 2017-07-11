# convert a TSV that was converted from a TraML file into a FASTA concatenated
# decoy database as Mayu input.

# Usage: ProtFDRTraML2Fasta.py library.tsv library.fasta

import sys, csv, argparse

def handle_args():
    usage = "" #usage: %prog --in \"files1 file2 file3 ...\" [options]" 
    usage += "\nThis program converts a TSV spectral library to FASTA"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", required=True, help = 'An OpenSWATH spectral library')
    parser.add_argument("--out", dest="outfile", required=True, help="Output file")
    parser.add_argument('--remove_nonunique', action='store_true', default=False, help="Remove non-unique ids (e.g. 2/id1/id2)")

    args = parser.parse_args(sys.argv[1:])
    return args

args = handle_args()

# infile = csv.reader(open(sys.argv[1]), delimiter="\t")
# outfile = open(sys.argv[2], "w")
infile = csv.reader(open(args.infile), delimiter="\t")
outfile = open(args.outfile, "w")

header = next(infile)
header_dict = dict([(h, i) for i,h in enumerate(header)])

pepseq_pos = header_dict["PeptideSequence"]
prot_pos = header_dict["ProteinName"]

class Protein(object):
    def __init__(self, name):
        self.peptides = set()
        self.name = name

    def add_peptide(self, peptide):
        self.peptides.update([peptide])

    def get_concat_peptides(self):
        return "".join(self.peptides)

protein_dic = {}
for line in infile:
    peptide = line[pepseq_pos]
    protein = line[prot_pos]
    if args.remove_nonunique and not (protein.startswith("DECOY_1/") or protein.startswith("1/")):
        continue
    if protein not in protein_dic:
        p = Protein(protein)
        protein_dic[protein] = p
    protein_dic[protein].add_peptide(peptide)


# output
for k in protein_dic:
    protein = protein_dic[k]
    outfile.write(">%s\n" % protein.name)
    outfile.write(protein.get_concat_peptides())
    outfile.write("\n")

