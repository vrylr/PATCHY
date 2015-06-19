#!/usr/bin/python

import sys

#helper script to generate list of positions to generate primers for


fwd_start = 1750
fwd_stop = 1806
fwd_increment = 3

rev_start = 1745
rev_stop = 1688
rev_increment = 3


fwd_outstring = ""
rev_outstring = ""

fwd_pos = fwd_start
rev_pos = rev_start

while( fwd_pos <= fwd_stop ):
	fwd_outstring = fwd_outstring + str( fwd_pos ) + " "
	fwd_pos = fwd_pos + fwd_increment

while( rev_pos >= rev_stop ):
        rev_outstring = rev_outstring + str( rev_pos ) + " "
        rev_pos = rev_pos - rev_increment

print fwd_outstring
print rev_outstring
