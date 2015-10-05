from __future__ import print_function
# Useage: python trafoXML_visualize.py input.trafoXML output_dir
import re, numpy, sys
from matplotlib.mlab import *
from matplotlib.pyplot import *

resdir = ''
file_in =  'small_002.trafoXML'
if len(sys.argv) > 1: file_in = sys.argv[1]
if len(sys.argv) > 2: resdir  = sys.argv[2]
f = open(file_in)
text = f.read()
f.close()

# parse the input file into pairs of x/y coordinates
pair_re = re.compile('<Pair from="([^ ]*)" to="([^ ]*)"/>')
x = []
y = []
for pair in pair_re.finditer(text):
    x.append( float(pair.group(1)) )
    y.append( float(pair.group(2)) )

# calculate least squares regression
A = np.vstack([x, np.ones(len(x))]).T
m, c = numpy.linalg.lstsq(A,y)[0]

print("Use linear fit", m, c)

# calculate resides
residues = []
predicted = []
for xx,yy in zip(x,y):
    residues.append( yy - (xx *m + c) )
    predicted.append( [xx, xx*m+c])

## plotting
plot(x,residues, 'ro' )
plot( [min(x), max(x)] , [0,0])
savefig(resdir+ 'trafo_residues_plot.pdf', format='pdf')
clf()

predicted.sort(key=lambda x: x[0], reverse=True)
plot(x,y, 'ro' )
plot( [xx[0] for xx in predicted], [yy[1] for yy in predicted] )
savefig(resdir+ 'rt_correlation.pdf', format='pdf')
clf()
