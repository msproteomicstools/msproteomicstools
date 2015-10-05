from msproteomicstoolslib.format.ProteinDB    import     ProteinDB  
import sys

def main(argv):
    
    fasta = argv[0]
    fasta_new = fasta + ".new"
    
    #Read fasta file
    protDB = ProteinDB()
    protDB.readFasta(fasta)
    
    for protein in protDB.proteinDictionary.values() :
        code1 = protein.code1
        code2 = protein.code2
        if len(code1) == 0 : code1 = code2 
        if len(code2) == 0 : code2 = code1
        protein.code1 = code2
        protein.code2 = code1

    protDB.writeFastaFile(fasta_new, chunksize = -1, format="sp")

if __name__ == '__main__' : 
    main(sys.argv[1:])