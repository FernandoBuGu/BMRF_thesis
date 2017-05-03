"""
Diego F. Montiel Gonzalez
P3 Excercise
"""
from Bio import Entrez
from sys import argv
import os.path
import p2

def readFileInput(filename):
    """
    Function: readFileInput
    Description: Get all the sequences of each identifier in the inputfile
    as a fasta format
    Input: Text file with the identifiers to look for: Ex. "NM_1799883"
    Output: List with all the sequences (similar to a fasta file)
    """
    Entrez.email = 'diego.montielgonzalez@wur.nl'
    identifiers = open(filename)
    gb_records = ''
    for iterator in identifiers.readline().rsplit():
        handle = Entrez.efetch(db="nucleotide",id=[iterator], rettype="gb")
        records = handle.read()
        gb_records = gb_records + records
    return gb_records

def generateGenBankFile(gb_records):
    """
    Function: generateGenBankFile
    to generate the genbankfile and then use the function of the p2 script
    to parse the genbank file
    Input: String with the Genbank records
    Output: Genbank file or an False boolean statement
    """
    genbank_file_name = "genbank_file.gb" 
    if not os.path.isfile(genbank_file_name):

        with open(genbank_file_name,"a+") as gb_file:
            gb_file.write(gb_records)
            if os.path.isfile(genbank_file_name):
                print "Succes!, Genbank File succesfully generated!"
                return genbank_file_name

            else:
                print "Failed to generated Genbank File" 
                return False
    else: 
        return genbank_file_name
        
if __name__ == '__main__':
    """
    Inputs: Script and text file with genbank identifiers
        Ex. python p3.py p3input.txt
    """
    sequence_array  = []
    filename        = argv[1] 
    gb_records      = readFileInput(filename)
    gb_file         = generateGenBankFile(gb_records)
    if  gb_file != False:
        accession, organism, sequence = p2.parseGenBankFile(gb_file)
        lenght_sequence = p2.getLenghtOfSequence(sequence)
        gc_content      = p2.getListOfGCContent(sequence)
        sequence_array  = p2.orderSequenceByGCContent(accession,organism,sequence,lenght_sequence,gc_content)     
        print p2.generateFastaFile(sequence_array)
        print p2.generateTextFile(sequence_array)
    else:
        print "Program failed...verified the outputs!"