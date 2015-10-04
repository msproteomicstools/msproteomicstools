import sys

f = open(sys.argv[1])
fout = open(sys.argv[2], "w")
header = next(f)
fout.write("preprocess_id\t" + header)
for i, line in enumerate(f):
    fout.write(str(i) + "\t" + line)
fout.close()
