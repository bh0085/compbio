#!/usr/bin/env python
import sys
import string
import random

#### INSTRUCTIONS FOR USE:
# call program as follows: ./gibbs.py <Motif Length> <Data File>
# make sure the gibbs.py is marked as executable: chmod +x gibbs.py

alphabet = ['A', 'G', 'C', 'T']

#### GibbsSampler:
#### 	INPUTS:	S - list of sequences
####		L - length of motif
####	OUTPUT:	PWM - 4xL list with frequencies of each base at each position
####                  Order of bases should be consistent with alphabet variable
def GibbsSampler(S,L):
    PWM = []
    for i in range(len(alphabet)):
        PWM.append([0]*L)

    ######### ADD YOUR CODE HERE ######

    ######### END OF YOUR CODE HERE #####
	
    return PWM

###### YOUR OWN FUNCTIONS HERE
# optional -- feel free to add your own functions if you want to


###### END OF YOUR FUNCTIONS

def main():
    L = int(sys.argv[1])
    datafile = sys.argv[2]
    S = readdata(datafile)
	
    P = GibbsSampler(S,L)
	
    print "    ", 
    for i in range(L):
        print "%-5d " % (i+1),
    print ""
	
    for j in range(len(alphabet)):
        print " %s " % alphabet[j], 
        for i in range(L):
            print " %5.3f" % P[j][i],
        print ""
	
def readdata(file):
    data = [];
    for line in open(file,'r'):
        data.append(line[0:-1])
    return data

main()

