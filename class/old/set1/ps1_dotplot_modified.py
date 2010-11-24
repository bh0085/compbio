#!/usr/bin/env python

#### 
# 6.047/6.878 - Problem Set 1 - string hashing/dotplots
#
# INSTRUCTIONS FOR USE:
# call program as follows: 
#  ./ps1-dotplot.py <FASTA 1> <FASTA 2> <PLOTFILE>
#     e.g. ./ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
#
# Make sure the ps1-dotplot.py is marked as executable: 
#     chmod +x ps1-dotplot.py
# or in windows with:       
#     python ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
# once you have put python in your path
#
#
# GNUPLOT 
# Gnuplot is used to generate plots for this program.  It is a common plotting
# program installed on most unix systems.  To use gnuplot on athena do the 
# following:
#
# athena% add gnu
#
# To test that it works do this:
#
# athena% gnuplot
# gnuplot> plot cos(x)
#
# you should then see a cosine plot appear.  
#
# Note: plotting.py and util.py must be in the same directory as this script.
# These files contain the code for generating plots.  You should not 
# worry about understanding any of the code contained within these files.
    # Much of it is copied from another project and is unrelated.
#



import sys, random
import plotting


def readSeq(filename):
    """reads in a FASTA sequence"""
    
    stream = open(filename)
    seq = []
    
    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())
    
    return "".join(seq)


def quality(hits):
    """determines the quality of a list of hits"""
    
    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000) 
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000
    
    goodhits = []
    
    for hit in hits:
        upper = slope1 * hit[0] + offset1
        lower = slope2 * hit[0] + offset2
        
        if lower < hit[1] < upper:
            goodhits.append(hit)
            
    return goodhits


def makeDotplot(filename, hits):
    """generate a dotplot from a list of hits
       filename may end in the following file extensions:
         *.ps, *.png, *.jpg
    """
    x, y = zip(* hits)
    
    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000) 
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    hits2 = quality(hits)
    print "%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits)))
    
    # create plot
    p = plotting.Gnuplot()
    p.enableOutput(False)
    p.plot(x, y, xlab="sequence 2", ylab="sequence 1")
    p.plotfunc(lambda x: slope1 * x + offset1, 0, 1e6, 1e5)
    p.plotfunc(lambda x: slope2 * x + offset2, 0, 1e6, 1e5)
            
    # set plot labels
    p.set(xmin=0, xmax=1e6, ymin=0, ymax=1e6)
    p.set(main="dotplot (%d hits, %.5f%% hits on diagonal)" %
          (len(hits), 100 * len(hits2) / float(len(hits))))
    p.enableOutput(True)
    
    # output plot
    p.save(filename)
    
    return p


inversion_map={'A':'T','T':'A','G':'C','C':'G'}
def inverted_appearance(key):

    #in order to increase the speed of this inversion, 
    #automatically invert all keys with more A than T
    #more C than G


    key_inv = ''
    

    for i in range(len(key)):
        v = key[-1 - i]
        key_inv += inversion_map[v]

    key_inv.replace('A','T')
    key_inv.replace('T','A')
    key_inv.replace('G','C')
    key_inv.replace('C','G')
    
    vals = {'G':3,'C':2,'T':1,'A':0}
    for i in range(len(key)):
        if vals[key_inv[i]] > vals[key[i]]:
            return key_inv
        elif vals[key_inv[i]] < vals[key[i]]:
            return key

    return key
    

def main(plotfile = 'dotplot.jpg',
         file1 = 'data/human-hoxa-region.fa', 
         file2 = 'data/mouse-hoxa-region.fa',
         skip = 1,
         delskip = 0,
         kmerlen = 30,
         allow_inversions = False):

    # NOTE to WINDOWS users:
    #   If you do not want to use the command line, comment out the command line
    #   parsing code and hard-code the input filenames.
    #
    # For example, use the following:
    # file1 = "human-hoxa-region.fa"
    # file2 = "mouse-hoxa-region.fa"
    # plotfile = "dotplot.jpg"

    if __name__ =='__main__':

        # parse command-line arguments
        if len(sys.argv) < 4:
            print "you must call program as:  "
            print "   ./ps1-2.py <FASTA 1> <FASTA 2> <PLOT FILE>" 
            print "   PLOT FILE may be *.ps, *.png, *.jpg"
            sys.exit(1) 
        file1 = sys.argv[1] 
        file2 = sys.argv[2]
        plotfile = sys.argv[3]
        
    
    # read sequences
    print "reading sequences"
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)
    
    
    
    # hash table for finding hits
    lookup = {}
    
    #MY CODE MODIFICATION FOR SKIPS etc
    #kmerlen: length of hash key (Including skips)
    #skip:    For simple skips, set a skip in key computation. (Input keyword 'skip')
    #delskip: For mismatch allowance, delete elements from the key with a skip. (Input keyword 'delskip')
    #allow_inversions:    Do we allow inversions by halving the hash table?


    # store sequence hashes in hash table
    print "hashing seq1..."
    for i in xrange(len(seq1) - kmerlen + 1):
        key = seq1[i:i+kmerlen:skip]
        if delskip > 1:
            key = list(key)
            del key[::delskip]
            key = ''.join(key)
        if allow_inversions:
            key = inverted_appearance(key)

        lookup.setdefault(key, []).append(i)
        if i % 10000 == 0 : print i/10000

                                   
    # look up hashes in hash table
    print "hashing seq2..."
    hits = []
    for i in xrange(len(seq2) - kmerlen + 1):
        key = seq2[i:i+kmerlen:skip]
        if delskip > 1:
            key = list(key)
            del key[::delskip]
            key = ''.join(key)
        if allow_inversions:
            key = inverted_appearance(key)

        # store hits to hits list
        for hit in lookup.get(key, []):
            hits.append((i, hit))
    
    #
    # hits should be a list of tuples
    # [(index1_in_seq2, index1_in_seq1),
    #  (index2_in_seq2, index2_in_seq1),
    #  ...]
    #
    
    print "%d hits found" % len(hits)    
    print "making plot..."
    p = makeDotplot(plotfile, hits)



if __name__ == '__main__':
    main()
