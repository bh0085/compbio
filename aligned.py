import re
import unafold
import numpy as np
import subprocess
import os

#Little routine to rip the dotted elts out of a sequence
def trimdots(str_in):
    l = np.array(list(str_in))
    str_out = ''.join(l[np.nonzero(l != '.')])
    return str_out

def consform(pattern):
    consform = ''
    for i in range(len(pattern)):
        if pattern[i] == 0:
            consform += '.'
        elif pattern[i] > i:
            consform += '<'
        else:
            consform += '>'
    return consform

#Print a phylip file in non-interleaved format
def print_phylip(aligned,fname):
    s = aligned.seqs
    ntaxa = len(s.keys())
    nalign= len(s.items()[0][1])
    f = open(fname,'w')
    f.write(' '.join([str(ntaxa),str(nalign)])+ '\n')
    for key, value in s.items():
        f.write('{0:30}'.format(key) + str(value).replace('.', '-')+'\n')
    f.close()
    
def make_tree(aligned):
    if not os.path.isdir('temp'):
        os.mkdir('temp')
    fname = 'temp/temp.phylip'
    print_phylip(aligned, fname)
    subprocess.call(["phyml -i " + fname],shell=True)
    
                
class Aligned:
    seqs ={}
    folds = {}
    consenus = None
    def __init__(self,alignments, folds = {}, consensus = None):
        self.seqs = alignments
        self.folds = folds
        self.consensus = consensus
        

    def getfold(self,key):
        seq = self.seqs[key]
        folded = unafold.fold(seq)
        self.folds[key] = folded

    def makefolds(self):
        keys = self.seqs.keys()
        for k in keys:
            if not self.folds.get(k):
                self.getfold(k)

    def printfolds(self):
        for k in self.seqs.keys():
            seq = self.seqs[k]
            fold = consform(self.folds[k][0])
            
            ctrim = np.array(list(self.consensus))
            ctrim = ''.join(ctrim[np.nonzero(np.array(list(seq)) != '.')])
            print fold

    def compsec(self,key):
        if not self.folds.get(key):
            self.getfold(key)
        fold = self.folds[key]
        pattern = fold[0]



        print 'Comparison of keyed sequence to consensus:'

        seq = self.seqs[key]
        ctrim = np.array(list(self.consensus))
        sl = np.array(list(seq))
        ctrim = ''.join(ctrim[np.nonzero(sl != '.')])
        

        consform = consform(pattern)
        print ctrim #trimdots(self.consensus)
        print consform #trimdots(consform)
            
        
        

def parse_stk():
    #Stockholm format files appear to consist
    #of a bunch of commented '#' headers followed
    #by groups of alignments seperated by doubleskips.

    #parse_stk advances throught he alignment groups,
    #adding sequences piecewise into a dictionary
    #keyed by species name.

    #Two additional keys sit in the dictionary
    #containing a secondary structure prediction
    #and something... else?

    fname = 'covary/5s.stk'
    seqdict  = {}
    extras = {}
    
    #Break the file up into groups of doubleskips
    grps = open(fname).read().split('\n\n')
    for g in grps:

        lines = g.split('\n')
        #Toss out the header.
        plus1 = lines[0]
        if plus1[0] == '#':
            continue

        #Special stuff for the tail group.
        nofs = -1
        while 1:
            if nofs < -1 * len(lines):
                print 'Something is suspicious here'
                raise Exception()
            elif not lines[nofs] or lines[nofs][0] != '#':
                nofs -= 1
            else:
                break
        #Kick out tail group junk.
        if nofs < -1:
            lines = lines[0:nofs + 1]
        last = lines.pop()
        seclast = lines.pop()
        
        for l in lines:
            #Fill a dictionary with sequences
            match = re.search(re.compile('([^\s]*)\s*([^\s]*)'),l)
            key = match.group(1)
            seqdict[key] = seqdict.get(key, '') +match.group(2)
            #And note where the second group of the match starts
            astart = match.start(2)
            
        lkey = last[0:astart].strip()
        lval = last[astart:].strip()
        skey = seclast[0:astart].strip()
        sval = seclast[astart:].strip()
        
        extras['cons'] = extras.get('cons','') + sval
        extras['rf'] = extras.get('rf','') + lval

    a = Aligned(seqdict)
    a.consensus = extras['cons'] 
    return a
        
    
    


def parse_fa():
    fname = 'covary/5s.stk'
    lines = open(fname).read()

    seq_re = re.compile('')
    raise Exception()
    
