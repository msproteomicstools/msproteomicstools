#!/usr/bin/env python
from __future__ import print_function
import sys, gzip, csv

if len(sys.argv) < 3:
    print("Usage: align2msstats.py input_feature_alignment.tsv[.gz] output_msstats_compatible.csv")
    sys.exit(1)
inarg = sys.argv[1]
outarg = sys.argv[2]

if inarg.endswith(".gz"):
    f = gzip.open(inarg, 'r')
else:
    f = open(inarg, 'r')
csvf = csv.DictReader(f, dialect=csv.excel_tab)

out = open(outarg, 'w')
out.write(
    "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Intensity,Condition,BioReplicate,Run\n")

for line in csvf:
    name = line["align_origfilename"].split("/")[-1]
    fas = line["aggr_Fragment_Annotation"].split(";")
    pas = line["aggr_Peak_Area"].split(";")
    for fa, pa in zip(fas, pas):
        out.write("%s,%s,%s,%s,NA,light,%s,%s,%s,%s\n" % (line["ProteinName"], line["FullPeptideName"], line["Charge"], fa, pa, line.get("Condition",""), line.get("BioReplicate",""), name))

f.close()
out.close()
