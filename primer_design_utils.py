#!/usr/bin/python

import sys
import re
import PrimerTMCalculator

def reverse_complement( sequence ):
    rev_comp_list = []
    length = len( sequence )
    for i in range( length ):
        cur_char = sequence[length - i -1]
        #print "at i=%s, char to check is %s"%(i, cur_char)
        if( cur_char == 'A' ):
            rev_comp_list.append('T')
        elif( cur_char == 'C' ):
            rev_comp_list.append('G')
        elif( cur_char == 'G' ):
            rev_comp_list.append('C')
        elif( cur_char == 'T' ):
            rev_comp_list.append('A')
        else:
            print "Error: When making reverse_complement, Sequence contained a '"+cur_char+"' character at position %s, quitting..."%(i+1)
            sys.exit()
    #print "rev complemnt: %s of length %s got turned into"%(sequence, length)
    #print rev_comp_list
    return "".join( rev_comp_list )


#function to generate a histogram from a list of values
#returns a string for now
def generate_histogram( list, increment ):
    returnstring = ''
    list.sort()
    low_val = list[0]
    cur_bin = (int( low_val / 1 ) ) * 1.0
    cur_bin_count = 0

    for i in list:
        #print "i is %s and cur_bin is %s"%(i, cur_bin)
        if (i - cur_bin ) <=increment:    #used to be if (i % cur_bin ) <=increment
            cur_bin_count = cur_bin_count + 1
        else:
            returnstring = returnstring + "%s-%s: %s\n"%(cur_bin, cur_bin + increment, cur_bin_count)
            cur_bin = float( cur_bin + increment )
            #print "maap setting curbin to %.2f" % cur_bin
            while (i - cur_bin ) >= increment:
                returnstring = returnstring + "%s-%s: %s\n"%(cur_bin, cur_bin + increment, 0)
                cur_bin = float( cur_bin + increment )
                #print "moop setting curbin to %.2f" % cur_bin
                
            cur_bin_count = 1
    returnstring = returnstring + "%s-%s: %s\n"%(cur_bin, cur_bin + increment, cur_bin_count)

    return returnstring



#very important function
#returns an integer that's supposed to
#reflect the quality of the primer
#return value of 0 means primer is 'perfect'
#for every problem, the return value gets implemented by 1
#current problems:
#1. TM stuff: superproblem TM too low, problem Tm too high
#2. at least 2 AT need to be in the last 5 nucleotides
#3. the last two nucleotides need to be a G or C
#4. superproblem: last nucleotide needs to be G or C, if not, increment by 10 and leave
#function could also be adapted to measure GC content and other stuff
def num_primer_problems( sequence, min_tm, max_noproblem_tm ):
    #print "checking problems of sequence %s" % sequence
    tm = PrimerTMCalculator.TM_Calculator().calculate_tm( sequence )

    if tm < min_tm:
        #print "tm is %.2f, returning 10" %tm
        return 10
    
    num_problems = 0
    if tm > max_noproblem_tm:
        num_problems = num_problems + 1

    
    length = len( sequence )
    last_5_nt = sequence[length - 5: length]
    #print "last_5_nt are %s" % last_5_nt
    if (last_5_nt[4] == 'A' ) or (last_5_nt[4] == 'T' ):
        #print "returning because last nt is a or t"
        return 10
    if (last_5_nt[3] == 'A' ) or (last_5_nt[3] == 'T' ): #in combination with previous if statement, ensures last 2 is g/c
        num_problems = num_problems + 1

    at_count = 0
    for i in range(4):
        if (last_5_nt[i] == 'A' ) or (last_5_nt[i] == 'T' ):
            at_count = at_count + 1
    if at_count < 2:
        num_problems = num_problems + 1
    #print "returning %s" % num_problems
    return num_problems


def read_insertion_site_list( filename, convert_list ):

    to_return = [];

    filename = filename.replace("\n","")
    fileh = open(filename, 'r')
    filelines = fileh.readlines()
    fileh.close()

    for line in filelines:
        line_items = line.split();
        for litem in line_items:
            this_pos = int( litem)
            if convert_list == 1:
                this_pos = 3 * this_pos + 382
            to_return.append( this_pos )
    return to_return


