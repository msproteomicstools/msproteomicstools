import pydot
import yaml
import sys

allowedTreeFormat = ["Mapped", "MappedFile", "MappedFileInput", "Raw"]
allowedDrawFormat = ["dot", "neato", "twopi", "circo", "fdp", "sfdp"]

if len(sys.argv) < 4:
    print "Usage: input.yaml output.pdf treeFormat drawFormat"
    print "  treeFormat is one of (%s)" % ", ".join(allowedTreeFormat)
    print "  drawFormat is one of (%s)" % ", ".join(allowedDrawFormat)
    print ""
    print "    this software takes a .yaml output file from a feature_alignment\n\
    run and draws the corresponding tree in a specific output format (pdf, svg, ...)"
    print "    requires the graphviz package to be installed and works on a linux machine"
    sys.exit()

tmpfile = "/tmp/tree_tmpfile.raw"
infile = sys.argv[1]
outfile = sys.argv[2]
treeFormat = sys.argv[3]
drawFormat = sys.argv[4]

if treeFormat not in allowedTreeFormat:
    print "Error, unrecognized format '%s'" % treeFormat

if drawFormat not in allowedDrawFormat:
    print "Error, unrecognized format '%s'" % treeFormat

stream = open(infile, "r")
docs = yaml.load(stream)
tree = docs["AlignedSwathRuns"]["Output"]["Tree"][treeFormat]
myFormat = outfile.split(".")[-1]

newtree = []
for a,b in tree:
    newtree.append( [a.split("/")[-1], b.split("/")[-1] ] )

tree = newtree

def get_graph(tree, labels, length, extra=""):
    res = ""
    extra = "overlap=scale"
    for i,edge in enumerate(tree):
        if labels is None or length is None:
            res += '"%s" -- "%s"\n' % ( edge[0], edge[1] )
        else:
            res += '"%s" -- "%s" [len=%s,label="%s"]\n' % ( edge[0], edge[1],  length[i], labels[i])
    return "graph G {\n%s\n%s\n}\n" % (res, extra)

gr = get_graph(tree, None, None)

with open(tmpfile, "w") as f:
    f.write(gr)

import subprocess
cmdline = "%s %s -T%s > %s" % (drawFormat, tmpfile, myFormat, outfile)
p = subprocess.Popen(
            cmdline,
            stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

# gr = pydot.graph_from_edges(tree)
# gr.write(outfile, prog=drawFormat, format=myFormat)

