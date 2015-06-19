#!/usr/bin/python

#note: converting serial cloner gb files to mac format
#in vi: %s/^M/\r/g
#see
#http://stackoverflow.com/questions/64749/m-character-at-end-of-lines
#http://stackoverflow.com/questions/71323/how-to-replace-a-character-for-a-newline-in-vim
#

import sys
import re

class GeneticCodeAAConverter:

    def __init__(self):

        self.aa2codons_dict = {
            'I' : ['ATT', 'ATC', 'ATA'],
            'L' : ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
            'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
            'F' : ['TTT', 'TTC'],
            'M' : ['ATG'],
            'C' : ['TGT', 'TGC'],
            'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
            'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
            'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
            'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
            'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'Y' : ['TAT', 'TAC'],
            'W' : ['TGG'],
            'Q' : ['CAA', 'CAG'],
            'N' : ['AAT', 'AAC'],
            'H' : ['CAT', 'CAC'],
            'E' : ['GAA', 'GAG'],
            'D' : ['GAT', 'GAC'],
            'K' : ['AAA', 'AAG'],
            'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            '*' : ['TAG', 'TAA' , 'TGA']
            }

        self.codon2aa_dict = {}
        for i in self.aa2codons_dict.keys():
            for codon in self.aa2codons_dict[ i ]:
                self.codon2aa_dict[ codon ] = i

    def aa2codons( self, aa ):
        if not aa2codons_dict.has_key( aa ):
            print "Error, no aa of 1-letter code "+aa+" known."
            sys.exit()
        return self.aa2codons_dict[aa]

    def codon2aa(self, codon ):
        if not self.codon2aa_dict.has_key( codon ):
            print "Error, no codon "+codon+" known."
            sys.exit()
        return self.codon2aa_dict[ codon ]


def read_dna_seq_from_gb_file( filename ):
    filename = filename.replace("\n","")
    fileh = open(filename, 'r')
    filelines = fileh.readlines()
    fileh.close()

    seq_read = ''
    seq_to_return = []

    sequence_reading_flag = 0
    for line in filelines:
        items = line.split()

        if( (len(items) == 1) and ( (items[0] == 'ORIGIN') or (items[0] == 'origin'))):
            sequence_reading_flag = 1
            continue
        if( (sequence_reading_flag == 1 ) and (items[0] == '//') ):
            sequence_reading_flag = 0
            continue

        if( sequence_reading_flag == 1 ):
            for i in range( len( items ) - 1 ):
                seq_read = seq_read + items[i+1] #+1 bc the first column is a number


    #sequence read in, now we do some quality control
    #and also convert everything to upper case letters
    #print "raw in is "
    #print seq_read
    for i in range( len(seq_read ) ):
        cur_char = seq_read[i]
        if( (cur_char == 'A') or (cur_char == 'C') or (cur_char == 'G') or (cur_char == 'T') ):
            seq_to_return.append( cur_char )
            continue
        
        if( cur_char == 'a' ):
            cur_char = 'A'
        elif( cur_char == 'c' ):
            cur_char = 'C'
        elif( cur_char == 'g' ):
            cur_char = 'G'
        elif( cur_char == 't' ):
            cur_char = 'T'
        else:
            print "Error: Sequence contained a '"+cur_char+"' character at position %s, quitting..."%(i+1)
            sys.exit()
        seq_to_return.append( cur_char )
        
    return "".join( seq_to_return )

'''
filename = ''

CommandArgs = sys.argv[1:]

for arg in CommandArgs:
    if arg == '-f':
        filename = CommandArgs[CommandArgs.index(arg)+1]

if( filename == '' ):
    print 'Error, please supply a file with option -f'
    sys.exit()



sequence = read_dna_seq_from_gb_file( filename )
print "the following sequence was obtained: "
print sequence
'''

