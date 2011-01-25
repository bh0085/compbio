#!/usr/bin/env python
import sys
import re

#In order to enumerate hexamers,
#convert numbers from 0 - 4095 to base 4.
def num2seq(num):
    seq = ''
    lookup = 'ATGC'
    
    current = num
    for dummy in range(6):
        base = current % 4
        seq += lookup[base]
        current/=4
    return seq

def parse_known_motifs(fname ='data/Yeast_known_motifs.txt'):
    data = open(fname).read()
    #omit the header line
    data = data[data.index('\n'):]
    #parse with a regular expression
    r = re.compile(r"(?P<name>\w+) (?P<sequence>\w+)",re.M)
    matches = list(re.finditer(r, data))
    parsed = {}
    parsed['names'] = map(lambda x: x.group('name'), matches)
    parsed['seqs'] = map(lambda x: x.group('sequence'), matches)
    return parsed


def main():
    if len(sys.argv) < 4:
        print "you must call program as:  "
        print "   ./hexamers.py <inter> <intercons> <output_file>"

    inter = open(sys.argv[1]).read()
    intercons = open(sys.argv[2]).read()
    bases = 'ATGC'
    
    #Create a hash and throw everything in.
    hexamers = {}
    for i in range(4096):
        hexamers[num2seq(i)] = [0,0]
    for i in range(len(inter) - 6):

        this_hex = inter[i: i + 6]
        #Check to make sure that we have a 'good' hexamer
        is_bad = False
        for h in this_hex:
            if h not in bases:
                is_bad = True
        if is_bad: continue

        #Check for conservation.
        this_cons = intercons[i:i+6]
        is_cons = True
        for item in this_cons:
            if item != '*':
                is_cons = False
        hexamers[this_hex][0] += 1
        if is_cons:
            hexamers[this_hex][1] += 1
 
    #Throw everything in to a list and save.
    hex_list = []
    for k, v in hexamers.iteritems():
        hex_list.append(  [k,v[0],v[1]])
    import pickle
    pickle.dump(hex_list,open(sys.argv[3],'w'))
       

def printout(data):

    #Calculate conservation fraction and append to raw data
    frac_cons = map(lambda x:float(x[2])/x[1], data)
    all_data = data
    for i in range(len(all_data)):
        all_data[i].append(frac_cons[i])

    #Compute sorted indices for each quantity of interest
    inds = range(len(all_data))
    freqs = map(lambda x: x[1], data)
    inds_freq = sorted(inds,key =  freqs.__getitem__, reverse = True)
    inds_frac = sorted(inds,key =  frac_cons.__getitem__, reverse = True)
    
    
    #DISPLAY CONSERVED REGIONS, ETC

    np = 50
    print 'setting np=' + str(np)
    #Print
    print 'Most frequently appearing sequences'
    print
    print '{0:20}{1:20}{2:20}{3:20}'.format('hexamer','appearances','conserved','frac')

    for j in range(np):
        i = inds_freq[j]
        print '{0:10}{1:20}{2:20}{3:20f}'.format(all_data[i][0],
                                                 all_data[i][1],
                                                 all_data[i][2],
                                                 all_data[i][3])
        
    print
    print 'Conservation frequency / appearance frequency'
    print '{0:20}{1:20}{2:20}{3:20}'.format('hexamer','appearances','conserved','frac')

    for j in range(np):
        i = inds_frac[j]
        print '{0:10}{1:20}{2:20}{3:20f}'.format(all_data[i][0],
                                                 all_data[i][1],
                                                 all_data[i][2],
                                                 all_data[i][3])
    


    #COMPUTE BINOMIAL STATISTICS
    #First, compute a gaussian approximation useful for large 
    #In order to determine the probability, r check the value of C/N
    #averaged at a few different N's (For high N, should converge to r)

    import numpy
    nt = len(all_data)
    nums = numpy.array(map( lambda x : x[1], all_data))
    cons = numpy.array(map( lambda x:  x[2], all_data))
    sortorder = numpy.argsort(nums)

    import math

    nbins = 10
    means = []
    for i in range(nbins):
        inds = sortorder[i * (nt / nbins) : (i+1)*(nt/nbins)]
        means.append((nums/cons)[inds].mean())

    #Output is:
    #    [10.209907509013952, 10.3382257012394, 10.284374773731084, 10.504155395495877, 
    #    10.123960612691468, 10.37140060977905, 10.141513101625067, 10.331578605316176, 10.365723361247474, 10.47246595465915]
    #I will use a value of p = 1/10.3 that seems to be close to the value towards upper end of the distribution in n.
        
    #The value of r at different appearance numbers is relatively consistent. For relatively large N, p values computed from
    #a gaussian distribution with occurence probably r should be an accurate representation.
    print means
    r = 1/(numpy.array(means).mean())
    pvals = []
    for i in range(nt):
        n = nums[i]
        c = cons[i]
        p = 0
        sig = math.sqrt(r * ( 1 -r) * n)
        nbar = r * n
        #compute the "integral"
        for x in range(c, n):
            p += (1.0 / (sig * math.sqrt(math.pi * 2))) * math.exp(- math.pow( float(x) - nbar,2) / (2 * math.pow(sig,2)))
        pvals.append(p)
                                   
    pval_sort = sorted(range(nt),key = pvals.__getitem__)
    
    print
    print 'Matches found:'
    print '{0:5} {4:8}{1:20}{2:10}{3:10} {5:10}{6:10} {7:10}'.format('n',
                                                                     'motif',
                                                                     'motif name',
                                                                     'match score',
                                                                     'hexamer',
                                                                     'number',
                                                                     'cons frac',
                                                                     'pval (%)',
                                                                     'bonf. 5%',
                                                                     'bonf. 20%')
    motifs = parse_known_motifs()
    motifs_found = []
    for j in range(np):
        i =  inds_frac[j]
        kmer = all_data[i][0]
        from compbio import sequtils
        
        num = all_data[i][1]
        frac = all_data[i][3]
        pval = pvals[i]

        found = sequtils.kmer_in_list(kmer, motifs['seqs'])
        if found[0] != -1 and found[1] > 3:
            print '{0:5} {4:8} {1:20} {2:10}{3:10} {5:10}{6:10} {7:10} {8:10}{9:10}'.format(j,
                                                                                      motifs['seqs'][found[0]],
                                                                                      motifs['names'][found[0]],
                                                                                      found[1],
                                                                                      kmer,
                                                                                      num,
                                                                                      frac,
                                                                                      pval * 100,
                                                                                      4096*pval < .05,
                                                                                      4096 * pval < .2)
        else:
            print '{0:5} {4:8} {1:20} {2:10}{3:10} {5:10}{6:10} {7:10} {8:10}{9:10}'.format(j,
                                                                                      'xxx','xxx','',
                                                                                      kmer,
                                                                                      num,
                                                                                      frac,
                                                                                      pval*100,
                                                                                      4096*pval < .05,
                                                                                      4096 * pval < .2)

    data_out = {'pvals':pvals,
                'hexamer':map(lambda x: x[0],all_data),
                'nums':map(lambda x: x[1],all_data),
                'cons':map(lambda x: x[2],all_data),
                'frac':map(lambda x: x[3],all_data)}
    return data_out
                
                               
            

        
if __name__ == "__main__":
    main()
