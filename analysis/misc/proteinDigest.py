from __future__ import print_function
import sys
import os
import csv
import getopt
import operator
import random

from msproteomicstoolslib.format.ProteinDB    import     ProteinDB  


def usage() :
    print("")
    print("proteinDBmasses.py")
    print ("-" * 15)
    print("This script retrieves the protein weights given a fasta file.")
    print("")
    print("Usage: ")
    print("python proteinDBmasses.py [options] fasta_file(s)")
    print("-h        Display this help")
    print("-e     enzyme    Enzyme used for in-silico digestion (peptide counting). Options: trypsin, Asp-N, Arg-C, Chymotrypsin, Lys-C, Lys-N. Default: trypsin")
    print("-l     pep-length    Minimum peptide length for the in-silico digestion. Default: 5")
    print("-m     missed-cleav  Allowed missed cleavages. Default: 0")
    print("")



def writeDigestFile(fastafile,enzyme, minPepLength = 5, missedCleavages = 0) :
    base = os.path.splitext(fastafile)[0]
    digest = base + '_digest.tsv'
    dynfile  = base + '_simulation.tsv'
    headers = ['protein_code','sequence']
    peptidelib = {} 
    
    db = ProteinDB()
    db.readFasta(fastafile)
    try :
        writer = csv.writer(open(digest,'w'), dialect='excel-tab')
    except :
        print("something went wrong while trying to write the file :" , massfile)
        sys.exit(1)
    
    writer.writerow(headers)
    protein_cnt = 0
    for code, protein in db.proteinDictionary.items() :
        protein_cnt += 1
        if protein_cnt % 5000 == 0 : print("%s proteins stored" % protein_cnt)
        peptides = protein.digest(enzyme, minLength = minPepLength, missedCleavages = missedCleavages)
        for peptide in peptides :
            wr_row = [protein.code2, peptide]
            writer.writerow(wr_row)
    
def main(argv) :

    trypsin = {'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : ['KP', 'RP']}
    Lys_N = {'terminus' : 'N' , 'cleave' : ['K'], 'exceptions' : []}
    Asp_N = {'terminus' : 'N' , 'cleave' : ['B','D'], 'exceptions' : []}
    Arg_C = {'terminus' : 'C' , 'cleave' : ['R'], 'exceptions' : ['RP']}
    Chymotrypsin = {'terminus' : 'C' , 'cleave' : ['F','Y','W','L'], 'exceptions' : ['FP','YP','WP','LP']}
    Lys_C = {'terminus' : 'C' , 'cleave' : ['K'], 'exceptions' : ['KP']}
    
    enzymes = {'trypsin' : trypsin , 'Lys-N' : Lys_N , 'Asp-N' : Asp_N, 'Arg-C' : Arg_C , 'Chymotrypsin' : Chymotrypsin , 'Lys_C' : Lys_C}
    enzyme = trypsin
    
    minPepLength = 5
    missedCleavages = 0
    
    #Get options
    try:
        opts, _ = getopt.getopt(argv, "he:l:m:",["help","enzyme","pep-length","missed-cleav"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    argsUsed = 0
    for opt,arg in opts:
        if opt in ("-h","--help") :
            usage()
            sys.exit()
        if opt in ("-e","--enzyme") :
            if arg not in enzymes :
                print("Error: Enzyme not recognized!") 
                print("Available enzyme options: " , [ key for key in enzymes.keys() ])
            enzyme = enzymes[arg]
            argsUsed += 2
        if opt in ("-l","pep-length") :
            minPepLength = int(arg)
            argsUsed += 2
        if opt in ("-m","missed-cleav") :
            missedCleavages = int(arg)
            argsUsed += 2

    fastafiles = argv[argsUsed:]
    for fastafile in fastafiles : 
        print("processing " , fastafile)
        #File exists?
        if not os.path.exists(fastafile) :
            print("This file: %s does not exist! It will be ignored." % fastafile)
            continue
        writeDigestFile(fastafile, enzyme, minPepLength, missedCleavages )
        
    pass

if __name__ == "__main__" :
    main(sys.argv[1:])