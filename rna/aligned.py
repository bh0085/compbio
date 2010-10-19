import re
import unafold
import numpy as np
import subprocess
import os

#Print a sequence output by unafold in the stockholm consensus format (<...>)
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

#Print a phylip file in non-interleaved phylip format
def print_phylip(aligned,fname):
    s = aligned.seqs
    ntaxa = len(s.keys())
    nalign= len(s.items()[0][1])
    f = open(fname,'w')
    f.write(' '.join([str(ntaxa),str(nalign)])+ '\n')
    for key, value in s.items():
        f.write('{0:30}'.format(key) + str(value).replace('.', '-')+'\n')
    f.close()

#Print an alignment to FASTA format
def print_fasta(aligned,fname):
    s = aligned.seqs
    f = open(fname,'w')
    for key, value in s.items():
        f.write('>'+key+'\n')
        f.write(str(value).replace('.', '-') + '\n')
    f.close()                
                
#Compute a tree for an alignment (currently uses phyml ...SLOW)
def make_tree(aligned):
    if not os.path.isdir('temp'):
        os.mkdir('temp')
    fname = 'temp/temp.phylip'
    print_phylip(aligned, fname)
    subprocess.call(["phyml -i " + fname],shell=True)
    
#Create a gblocked alignment from a gapped one.
#Good for creating parsimonious phylo trees out of noisy data.
def gblocked_alignment(aligned):
    if not os.path.isdir('temp'):
        os.mkdir('temp')
    fname = 'temp/temp.fa'
    print_fasta(aligned, fname)
    os.chdir('temp')

    #gblocks command line:
    #-t    sequency type d = dna
    #-b1   minimum number of sequence for conserved
    #-b2   minimum number of sequence for flanking
    #-b3   maximum number of contig. nonconserved
    #-b4   minimum block length
    #-b5   allowed gaps  n (none), h (with half??), a (all)
    subprocess.call('Gblocks temp.fa -t=d -b3=8 -b4=4 -b5=h -v=250',
                    shell = True)
    print 'gblocking with default params: '
    print 'b1 = 50%+1'
    print 'b2 = 85%'
    print 'b3 = 8'
    print 'b4 = 4'
    print 'b5 = h'

    gba = parse_fa('temp.fa-gb') 
    os.chdir('..')
    return gba
   
#Geneeric class for holding aligned RNA/DNA and maybe some secondary
#structure.
class Aligned:
    seqs ={}
    folds = {}
    consenus = None

    #Init with an alignment etc.
    def __init__(self,alignments, folds = {}, consensus = None):
        self.seqs = alignments
        self.folds = folds
        self.consensus = consensus
        

    #Generate folded secstruct for an element in the alignment.
    def getfold(self,key):
        seq = self.seqs[key]
        folded = unafold.fold(seq)
        self.folds[key] = folded

    #Generate folded secstructs for each element in the alignment.
    def makefolds(self):
        keys = self.seqs.keys()
        for k in keys:
            if not self.folds.get(k):
                self.getfold(k)

    #Print a list fo all of the folded sequences
    def printfolds(self):
        for k in self.seqs.keys():
            seq = self.seqs[k]
            fold = consform(self.folds[k][0])
            
            ctrim = np.array(list(self.consensus))
            ctrim = ''.join(ctrim[np.nonzero(np.array(list(seq)) != '.')])
            print fold

    #Compare a sequence in the alignemnt to the consensus
    def compsec(self,key):
        if not self.folds.get(key):
            self.getfold(key)
        fold = self.folds[key]
        pattern = fold[0]
        print 'Comparison of keyed sequence to consensus:'
        seq = self.seqs[key]
        
        #Some juju to put the folded sequence and the
        #consensus sequence in the same format for viewing
        ctrim = np.array(list(self.consensus))
        sl = np.array(list(seq))
        ctrim = ''.join(ctrim[np.nonzero(sl != '.')])
        consform = consform(pattern)

        
        print ctrim #trimdots(self.consensus)
        print consform #trimdots(consform)
        
#Parse a stockholm file (e.g an rfam alignment of families
#Return an 'Aligned' object.
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

    #Create an aligned object
    a = Aligned(seqdict)
    a.consensus = extras['cons'] 
    return a
        
    
    

#Parses a fasta file and returns an 'Aligned' object
def parse_fa(fname):
    lines = open(fname).read().split('\n')

    #Create a dictionary with name/sequence pairs.
    seqdict = {}
    current_key = None
    for l in lines:
        if not l:
            continue
        elif not current_key and l[0] != '>':
            raise Exception()
        elif l[0] == '>':
            current_key = l[1:].strip()
        else:
            seqdict[current_key] = seqdict.get(current_key,'') + l.strip().replace(' ','')

    #Create an aligned object
    a = Aligned(seqdict)
    return a

    
