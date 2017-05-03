#!usr/bin/env python

"""
Names: Guus van de Steeg, Diego Montiel, Oscar Llorian
Student: ####
"""

from sys import argv
from numpy import mean

def parse_fastq(filename):
    """
    Function that returns each entry as a list with each sequence and its report

    Input: fastq file    
    """
    container = []
    fastq = open(filename)
    x = 0
    for line in fastq:        
        if x == 4:
            yield container
            print container
            container = []
            x = 0
        line = line.rstrip()
        if x != 2:
            container += [line]
        x += 1    
    fastq.close()

def fastq_translator(record):
    """Returns a string of numbers in ASCII format
    
    Input: List of 3 elements from parse_fastq()
    """
    ascii = ''
    quality = record[2]
    for element in quality:
        ascii += str((ord(element)-64)) + ' '
    ascii = ascii[:-1]
    record[2] = ascii
    return record
    
def length(yielded_record):
    """Returns the length of the sequence from the fastq record
    
    Input List of 3 elements from parse_fastq()
    """
    lengths = []
    for record in yielded_record:
        lengths += [len(record[1])]
    return min(lengths), mean(lengths), max(lengths)
        

if __name__ == '__main__':
    #test = parse_fastq(argv[1])
    test = parse_fastq('tomatosample.fq')
    for i in test:
        i = fastq_translator(i)
        break
    print length(test)
        