import csv
import xlwt


def getwriter(matrix_outfile):
    if matrix_outfile.endswith("xls"):
        matrix_writer = XlsWriter(matrix_outfile)
    elif matrix_outfile.endswith("tsv"):
        matrix_writer = CsvWriter(matrix_outfile,delim="\t")
    elif matrix_outfile.endswith("csv"):
        matrix_writer = CsvWriter(matrix_outfile, delim=",")
    else:
        raise Exception("Unknown matrix extension, must be .xls or .tsv")

    return matrix_writer


class CsvWriter():
    def __init__(self, outfile,delim="\t"):
        self.outfile = outfile
        self.f = open(outfile, "w")
        self.tsv = csv.writer(self.f, delimiter=delim)
        self.line = []

    def write(self, entry, color="ignored"):
        self.line.append(entry)

    def newline(self):
        self.tsv.writerow(self.line)
        self.line = []

    def __del__(self):
        print "Closing ", self.outfile
        self.f.close()


class XlsWriter():
    def __init__(self, outfile):
        self.outfile = outfile
        self.workbook = xlwt.Workbook()
        self.worksheet = self.workbook.add_sheet('0')
        self.col = 0
        self.row = 0
        self.colors = {'r': xlwt.easyxf('font: color red;'), 'b': xlwt.easyxf('font: color blue;'), 'd': xlwt.easyxf()}

    def write(self, entry, color='d'):
        thestyle = self.colors[color]
        self.worksheet.write(self.row, self.col, entry, thestyle)
        self.col += 1

    def newline(self):
        self.row += 1
        self.col = 0

    def __del__(self):
        print "Writing out", self.outfile
        self.workbook.save(self.outfile)