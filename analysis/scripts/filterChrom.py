import sys
sys.path.append("/media/data/var/openms/gitsvn/OpenMS/build4/pyOpenMS")
import pyopenms

if len(sys.argv) < 4:
  print "A small program that will filter chromatogram mzML file by their native id" 
  print "Usage: python filterChrom.py infile.mzML outfile.mzML chrom_id" 
  exit()

infile = sys.argv[1]
outfile = sys.argv[2]
filter_criteria = sys.argv[3]

print "Will filter with criteria ", filter_criteria

exp = pyopenms.MSExperiment()
pyopenms.MzMLFile().load(infile, exp)
exp2 = exp
exp2.clear(False)
chroms = exp2.getChromatograms()
chroms_out = []
import re
for c in chroms:
    #if c.getNativeID().find(filter_criteria) != -1:
    if re.search(filter_criteria, c.getNativeID()):
        chroms_out.append(c)

chroms_out.sort(lambda x,y: cmp(x.getNativeID(), y.getNativeID()))

exp2.setChromatograms(chroms_out)
pyopenms.MzMLFile().store(outfile, exp2)


