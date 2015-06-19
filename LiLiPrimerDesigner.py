#!/usr/bin/python

#script to read in a DNA sequence, list of fwd and rev positions,
#and then design primers that would amplify the plasmid at these positions
#for later blunt-end self-circularization
#Florian Richter, Moeglich Lab, oct 2014, florian.richter.1@hu-berlin.de,

import sys
import re

from ReadDNASeqs import read_dna_seq_from_gb_file
import PrimerTMCalculator
from primer_design_utils import reverse_complement, generate_histogram,num_primer_problems, read_insertion_site_list


MIN_Tm = 55.0 #the minimum annealing temp for designed primers
MAX_NON_PROBLEMATIC_TM = 69
MAX_P_LENGTH = 30 #maximum allowed length of the primer
MIN_P_LENGTH = 16
DESIRED_TM = 60.0 #desired 
TM_DIFF_PROBLEMATIC_INCREMENT = 2.5

ADD_LINKERS = 0
FWD_LINKER = "CCGGTACTGGCGGGTCGACT"
REV_LINKER = "GATCCACCTGCTGCTGCTGC"

#workhose function, scans through the sequence
#at site from min to max length, records problems
#for each, and returns the one with least problems
#orientation +1 = fwd, -1 = rev
#returns a sequence and a number of problems
def design_lili_primer( sequence, site, orientation ):
    options = []
    length = MIN_P_LENGTH

    best_p = ["",100, 0]

    while length <= MAX_P_LENGTH:
        candidate = ""
        if orientation == 1:
            candidate = sequence[site -1 : site+length -1 ]
        elif orientation == -1:
            candidate = reverse_complement( sequence[site -1 -length:site -1] )

        num_raw_problems = num_primer_problems( candidate, MIN_Tm, MAX_NON_PROBLEMATIC_TM )
        this_tm = PrimerTMCalculator.TM_Calculator().calculate_tm( candidate )

        tm_dev_problems = int( abs( DESIRED_TM - this_tm ) / TM_DIFF_PROBLEMATIC_INCREMENT )

        #print "Candidate %s for site %s has tm of %.2f, %s raw and %s tm_dev_problems" % (candidate, site, this_tm, num_raw_problems, tm_dev_problems)
        tot_problems = num_raw_problems  + tm_dev_problems
        
        if tot_problems < best_p[1]:
            best_p[0] = candidate
            best_p[1] = tot_problems
            best_p[2] = this_tm

        if tot_problems == 0:
            break
        length = length + 1

        #print "returning %s (length %s) for site %s with tm of %.2f and %s total problems" % (best_p[0], len(best_p[0]),  site, best_p[2] , best_p[1] )
    return best_p


#we're interested in average lengths and TM diff distribution
def calculate_statistics( primers ):
    lengths = []
    melt_temps = []

    for primer in primers:
        lengths.append( len( primer[0] ) )
        melt_temps.append( primer[2] )

    len_hist = generate_histogram( lengths, 1.0)
    tm_hist = generate_histogram( melt_temps, 2.5)

    return "Length distribution:\n" + len_hist + "\nTm distribution:\n" + tm_hist



seq_filename = ''
fwd_fname =''
rev_fname =''
fwd_list = []
rev_list = []

CommandArgs = sys.argv[1:]

for arg in CommandArgs:
    if arg == '-f':
        seq_filename = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-fwd_pos':
        fwd_fname =  CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-rev_pos':
        rev_fname = CommandArgs[CommandArgs.index(arg)+1] 
    elif arg == '-add_linkers':
        ADD_LINKERS = 1
    elif arg == '-desired_tm':
        DESIRED_TM = float( CommandArgs[CommandArgs.index(arg)+1] )

if( (seq_filename == '') or (fwd_fname == '') or (rev_fname == '') ):
     print 'Error, please specify sequence filename and fwd_pos and rev_pos filenames'
     sys.exit()


fwd_list = read_insertion_site_list( fwd_fname, 0 )
rev_list = read_insertion_site_list( rev_fname, 0 )

sequence = read_dna_seq_from_gb_file( seq_filename )
#print "obtained sequence:"
#print sequence

fwd_primers = []
rev_primers = []

#design primers
for fpos in fwd_list:
    fwd_primers.append( design_lili_primer( sequence, fpos, 1 ) )

for rpos in rev_list:
    rev_primers.append( design_lili_primer( sequence, rpos, -1 ) )

#now let's do some statistics
stat_string = calculate_statistics( fwd_primers + rev_primers )

outstring = "\n" + stat_string
outstring = outstring + "\nfwd_primers:\n"

if ADD_LINKERS == 1:
    for fprimer in fwd_primers:
        fprimer[0] = FWD_LINKER + fprimer[0]
    for rprimer in rev_primers:
        rprimer[0] = REV_LINKER + rprimer[0]

for i in range (len( fwd_primers )):
    outstring = outstring + "flib1_%s_fwd: %s, len=%s, tm=%.2f, num_probs=%s \n" % (fwd_list[i], fwd_primers[i][0], len(fwd_primers[i][0]), fwd_primers[i][2], fwd_primers[i][1] )

outstring = outstring + "\nrev_primers:\n"
for i in range (len( rev_primers )):
    outstring = outstring + "flib1_%s_rev: %s, len=%s, tm=%.2f, num_probs=%s \n" % (rev_list[i], rev_primers[i][0], len(rev_primers[i][0]), rev_primers[i][2], rev_primers[i][1] )
    
print outstring

    
