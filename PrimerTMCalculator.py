#!/usr/bin/python


#class to calculate primer TM based on breslauer pnas 83, 3746-3750 (1986)
#and http://www.basic.northwestern.edu/biotools/OligoCalc.html

import sys
import re
import math


class TM_Calculator:

    def __init__(self):
        #first we build table2 from the breslauer paper
        #slightly different implementation, we make an
        #entry for all 16 combinations of two bases
        self.twobase_data = {}
        #1. is deltaH in kcal/mol, 2. is deltaS in cal/K mol
        self.twobase_data['AA'] = [9.1, 24.0]
        self.twobase_data['AC'] = [6.5, 17.3]
        self.twobase_data['AG'] = [7.8, 20.8]
        self.twobase_data['AT'] = [8.6, 23.9]
        
        self.twobase_data['CA'] = [5.8, 12.9]
        self.twobase_data['CC'] = [11.0, 26.6]
        self.twobase_data['CG'] = [11.9, 27.8]
        self.twobase_data['CT'] = [7.8, 20.8]

        self.twobase_data['GA'] = [5.6, 13.5]
        self.twobase_data['GC'] = [11.1, 26.7]
        self.twobase_data['GG'] = [11.0, 26.6]
        self.twobase_data['GT'] = [6.5, 17.3]

        self.twobase_data['TA'] = [6.0, 16.9]
        self.twobase_data['TC'] = [5.6, 13.5]
        self.twobase_data['TG'] = [5.8, 12.9]
        self.twobase_data['TT'] = [9.1, 24.0]

        self.temp = 298
        self.initiation_const = -4.7 #-5.4 according to paper, -4.7 empirically more in line w/ thermo phusion website
        self.primer_conc = 0.0000005 #0.5uM
        self.salt_conc = 0.05 #50mM
        self.R = 0.0019872  #kcal /K mol
        #print "Enthalpy/Entropy table initiated."
        #print self.twobase_data


    def calculate_tm( self, sequence):
        sum_dH = 0
        sum_dS = 0
        seq_len = len( sequence ) - 1
        for i in range( seq_len ):
            pair=sequence[i:i+2]
            #print "quering "+pair+"..."
            if( self.twobase_data.has_key( pair ) ):
                pair_data = self.twobase_data[pair]
                sum_dH = sum_dH + pair_data[0]
                sum_dS = sum_dS + pair_data[1]
            else:
                print "Error: pair "+pair+" doesn't have data in the table."
                sys.exit()
        sum_dS = sum_dS / 1000
        dG = self.initiation_const + sum_dH - (self.temp * sum_dS)

        melt_temp = ( (sum_dH + self.initiation_const) / (sum_dS + self.R * math.log(1/self.primer_conc)) ) + (16.6 * math.log10(self.salt_conc)) - 273.15
        
        #print "sequence "+sequence+" has dH of "+str(sum_dH)+" and dS of %.3f"%sum_dS +" and a Tm of %.3f"%melt_temp + " and a deltaG of %.2f"%dG
        return melt_temp
            


sequence = ''
'''
CommandArgs = sys.argv[1:]

for arg in CommandArgs:
    if arg == '-s':
        sequence = CommandArgs[CommandArgs.index(arg)+1]

if( sequence == '' ):
     print 'Error, please supply a dna seq with option -s'
     sys.exit()

     #print "instantiating TM_calculator... "
calculator = TM_Calculator()

#print "done.."

calculator.calculate_tm( sequence )

'''

'''

#breslauer data
        self.twobase_data['AA'] = [9.1, 24.0]
        self.twobase_data['AC'] = [6.5, 17.3]
        self.twobase_data['AG'] = [7.8, 20.8]
        self.twobase_data['AT'] = [8.6, 23.9]
        
        self.twobase_data['CA'] = [5.8, 12.9]
        self.twobase_data['CC'] = [11.0, 26.6]
        self.twobase_data['CG'] = [11.9, 27.8]
        self.twobase_data['CT'] = [7.8, 20.8]

        self.twobase_data['GA'] = [5.6, 13.5]
        self.twobase_data['GC'] = [11.1, 26.7]
        self.twobase_data['GG'] = [11.0, 26.6]
        self.twobase_data['GT'] = [6.5, 17.3]

        self.twobase_data['TA'] = [6.0, 16.9]
        self.twobase_data['TC'] = [5.6, 13.5]
        self.twobase_data['TG'] = [5.8, 12.9]
        self.twobase_data['TT'] = [9.1, 24.0]

#sugimoto data
        self.twobase_data['AA'] = [8.0, 21.9]
        self.twobase_data['AC'] = [9.4, 25.5]
        self.twobase_data['AG'] = [6.6, 16.4]
        self.twobase_data['AT'] = [5.6, 15.2]
        
        self.twobase_data['CA'] = [8.2, 21.0]
        self.twobase_data['CC'] = [10.9, 28.4]
        self.twobase_data['CG'] = [11.8, 29.0]
        self.twobase_data['CT'] = [6.6, 16.4]

        self.twobase_data['GA'] = [8.8, 23.5]
        self.twobase_data['GC'] = [10.5, 26.4]
        self.twobase_data['GG'] = [10.9, 28.4]
        self.twobase_data['GT'] = [9.4, 25.5]

        self.twobase_data['TA'] = [6.6, 18.4]
        self.twobase_data['TC'] = [8.8, 23.5]
        self.twobase_data['TG'] = [8.2, 21.0]
        self.twobase_data['TT'] = [8.0, 21.9]

'''
