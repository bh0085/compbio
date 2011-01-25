import re
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import utils.seismic as seismic
import utils.plots as myplots

import utils.colors as mycolors

def parse():
    r =Reads('data/yeast_reads_aligned.txt')
    return r

class Reads():
    def __init__(self,fname):
        reads = open(fname).read()
        r = re.compile('^(\S+)\s+(\S+)\s+(\S+)',re.M)
        matches = re.finditer(r, reads)
        
        self.read_seqs = []
        self.read_refs = []
        self.read_starts = []
        self.ref_len = []
        self.read_len = []

        for m in matches:
            self.read_seqs.append(m.group(1))
            self.read_refs.append(m.group(2))
            self.read_starts.append(int(m.group(3)))
            self.read_len.append(len(m.group(1)))
            self.ref_len.append(len(m.group(2)))

        self.coverage = get_coverage(self)
        self.ref = get_ref(self)

        self.lens = array(self.ref_len,int)
        self.starts = array(self.read_starts,int)
        self.finishes = self.lens + self.starts

        
        self.read_idx_dict, self.reads = get_reads(self)
        self.refarray = get_refarray(self)
        self.error = q4b(self)

def show_reads(r):
    reflen = array(r.ref_len)
    readlen = array(r.read_len)
    #good, the numbers are always correct!
    print len(nonzero(not_equal(reflen,readlen))[0])
    
def get_coverage(r):
    """plot the coverage across the genome"""
    lens = array(r.ref_len)
    starts = array(r.read_starts)
    finishes = lens + starts
    m = max(finishes)
    coverage = zeros(m,float)

    n_reads = len(lens)

    for i in range(n_reads):
        coverage[starts[i]:finishes[i]] += 1

    #f = plt.figure(0)
    #f.clear()
    #seismic.seismic(f,[coverage[0:1000]],subplot = '111')
    
    return coverage

def get_ref(r):
    """Compile the reference genome from individual reference reads"""
    lens = array(r.ref_len,int)
    starts = array(r.read_starts,int)
    finishes = lens + starts
    n = len(starts)
    
    s0 = min(starts)
    f0 = max(finishes)
    
    ref = ['?' for i in range(f0)]
    for i in range(n):
        ref[starts[i]:finishes[i]] = r.read_refs[i]

    return ref

def get_refarray(r):
    d = r.read_idx_dict
    ref = r.ref
    refarray = zeros((len(ref),len(d.keys())),int)
    d2 = dict(d)
    d2['?'] = 4
    for i in range(len(ref)):
        refarray[i][d2[ref[i]]] = 1
    return refarray

def get_reads(r):
    """Get coverage counts for each letter in the reads"""
    idx_dict = {'A':0,'T':1,'C':2,'G':3,'N':4}
    n = len(r.starts)
    ni = len(idx_dict.keys())

    reads = zeros((max(r.finishes),ni),int)
    for i in range(n):
        start = r.starts[i]
        for j in range(r.lens[i]):
            reads[start+j][idx_dict[r.read_seqs[i][j]]] += 1

    return idx_dict, reads
    
def inds_at_which(r, which='mismatch'):
    if which == 'mismatch':
        inds = nonzero( np.sum(( 1-r.refarray) * r.reads, 1))[0]
    elif which == 'highcoverage':
        inds = nonzero( greater( r.coverage, 10))[0]
    else:
        raise Exception( " which type not yet handled")
    return inds

def show_reads(r,inds, fig = None, subplot = 111, flagged = None):
    ni = shape(r.reads)[1]
    ct = mycolors.getct(ni)
    if not fig: 
        fig = plt.figure(0)
        fig.clear()

    ax = fig.add_subplot(subplot)
    seismic.seismic(fig, r.reads[inds].T,
                   ax = ax, colors = ct,continuous = False,
                    labels = r.read_idx_dict.keys(),yunits = 'reads'
                    ,xunits = 'position from '+str(inds[0]))
    
    refarray = r.refarray
    seismic.seismic(fig, r.reads[inds].T * refarray[inds].T,
                    ax = ax, colors = 'w'*ni,continuous = False,
                    label_x = False, label_y = False)


    if flagged != None:
        seismic.seismic(fig, r.reads[inds].T * flagged,
                        ax = ax, colors = ['orange']* 5,continuous = False,
                        label_x = False, label_y = False)       
    
    return ax

#show a region of high mismatch

def get_seq_counts(r):
    inds = range(len(r.reads))
    maxidxs = []
    
    snps = zeros(len(inds)) -1
    ncorrect = zeros(len(inds))
    nwrong = zeros(len(inds))

    eps = .005

    d = r.read_idx_dict
    d2 = dict(d)
    d2['?'] = 4

    for i in inds:
        maxidx = argmax(r.reads[i])
        if maxidx == d2[r.ref[i]]:
            continue
        ncorrect[i] =  r.reads[i,maxidx]
        nwrong[i] = r.coverage[i] - ncorrect[i]

    min_log_confidence = log(238000 * 20)
    logprob = ( nwrong - ncorrect) * log(eps)
    wherecertain = nonzero(greater(logprob,min_log_confidence))[0]


    return logprob, wherecertain, maxidx

def read_position_errors(r):
    refs = r.ref

    maxlen = max(r.read_len)
    start_hash = zeros(maxlen)
    start_hash_t = zeros(maxlen)
    finish_hash = zeros(maxlen)
    finish_hash_t = zeros(maxlen)
    
    for i  in range(len(r.read_seqs)):
        rs = r.read_seqs[i]
        ofs = r.read_starts[i]
        l = len(rs)
        for j  in range(l):
            from_start = j
            from_finish = l - 1- j
            
            finish_hash_t[from_finish] += 1
            start_hash_t[from_start] += 1
            if rs[j] != refs[ofs + j]:
                finish_hash[from_finish] += 1
                start_hash[from_start] += 1
        
    
    print start_hash, start_hash_t
    return start_hash, start_hash_t, finish_hash, finish_hash_t

def show_pos_errors(h, ht):
    f = plt.figure(0)
    f.clear()
    ax = seismic.seismic(f,[h/ht],
                         y_marked = .02,
                       xunits = 'read position from start',
                       yunits = 'fraction')
    myplots.maketitle(ax, 'Read Errors by Position', 'with fractional errors on y axis')

    print mean(h/ht)

def show_certainties(r, lp, wc,maxi):
    maxi =argmax(r.reads,1)
    idxbases = [[]]*5
    for k,v in r.read_idx_dict.iteritems():
        idxbases[v] = k

    f = plt.figure(0)
    f.clear()
    ax = f.add_subplot(111)

    srtidx = argsort(lp)[::-1]
    print r.coverage[srtidx[0:1000]]
    ax.plot(lp[srtidx[0:1000]])
    
    print 'ref base & snp base & hits & misses & prob & context  \\\\'
    for i in range(50):
        idx = srtidx[i]
        prob =lp[srtidx[i]]
        rextr = r.ref[idx-5:idx+5]

        maxbase = idxbases[maxi[idx]]
        oldbase = r.ref[idx]
        hits = r.reads[idx][maxi[idx]]
        misses = r.coverage[idx]-hits

        print ' & '.join(["%s" % elt for elt in [maxbase,oldbase,hits,misses,round(prob,3),''.join(rextr)]]) +'''\\\\'''
    
    

def check_motifs():
    s = 'AAACACTTGGTGATTTGAATGTAATTTGAACGTTTAAAAAATTCCAAGGGAATATTACTGTTTCGGGAATATAACGTTTGATCGCTAGCGCACCCAGATGAGAATCCGCATGGATGGGTTGTTAGTATTTCTATTGTGAAACCGCTTTGTTCTGGAGGGCGAGAAAAAGTTACGGTCACTCTATCTTTTTTATAATTTGGTCT'
    import hexamers
    motifs = hexamers.parse_known_motifs('/Users/bh0085/Programming/Python/compbio/data/Yeast_known_motifs.txt')
    import utils.sequtils as seq
    print s
    m = seq.motifs_in_seq(s,motifs)
    
    

def show_zoom(r,matches = True):
    """Produce the plot 4c showing a small (presumably well linked) region where the sequence differs from the reference"""
    mm_inds = inds_at_which(r, 'mismatch')
    mmi0 = 1487
    mmi1 = 1525
    i0 = mm_inds[mmi0]
    i1 = mm_inds[mmi1]
    inds_all = range(i0, i1)
    
    

    ref = [r.ref[i] for i in inds_all]
    print 'refseq', ''.join(ref)
    refseq = ''.join(ref)

    motif = 'GATGAG'
    print 'MOTIF: ', refseq.index(motif)

    minds = arange(refseq.index(motif),refseq.index(motif)+7)
    gtginds = inds_all
    gtgon = array(gtginds) * 0.0
    gtgon[minds] = 1

    #gtginds = inds_all[refseq.index(motif):refseq.index(motif)+7]

    f0 = plt.figure(0)
    f0.clear()
    ax0 = show_reads(r,  inds_all, fig = f0, subplot = '311')
    ax1 = show_reads(r,  mm_inds[mmi0:mmi1],fig = f0,subplot ='312')
    ax2 = show_reads(r,  gtginds,fig = f0,subplot ='313', flagged = gtgon)

    myplots.color_legend(f0, mycolors.getct(5), r.read_idx_dict.keys())
    myplots.maketitle(ax0, 'Region under positive selection', 'Bars are histograms of reads at a position, colored bars are reads mismatched to reference position. \n("N" mismatches are ambiguous reads)')
    myplots.maketitle(ax1  ,'', 'Same region, only reference mismatches shown')
    myplots.maketitle(ax2  ,'', 'Same region, yeast motif GATGAG highlighted in orange')
    
def get_biases(r):
    refcounts = {'A':0,
                 'T':0,
                 'C':0,
                 'G':0}

    readcounts = {'A':0,
                 'T':0,
                 'C':0,
                 'G':0}

    for rs in read_seqs:
        for k in readcounts.keys():
            readcounts[k] += rs.count(k)
    for rs in read_refs:
        for k in radcounts.keys():
            refcounts[k] += rs.count(k)

    lens = array(r.ref_len)
    starts = array(r.read_starts)
    finishes = lens + starts
    
    total_ref = max(finishes) - min(starts) + 1
    total_reads = max
    


def q4b(r):
    high = inds_at_which(r, 'highcoverage')
    read_arr = r.reads[high,:]
    
    all_cons = []
    all_rights = []
    all_wrongs = []
    error_matrix = zeros((5,5))
    for i in range(len(read_arr)):
        cons = np.argmax(read_arr[i])
        all_cons.append(cons)
        other_idx = range(5); other_idx.remove(cons)
        all_wrongs.append(sum(read_arr[i,other_idx]))
        all_rights.append(read_arr[i,cons])
        error_matrix[cons,:] += read_arr[i]
    
    
    all_rights = array(all_rights)
    all_wrongs = array(all_wrongs)
    all_cons = array(all_cons)

    print 'Error rate over all high coverage areas'
    print float(np.sum(all_wrongs))/float(np.sum(all_rights + all_wrongs))
    print 'Error transition matrix'
    print error_matrix
    print 'Error probabilities'
    em2 = array(error_matrix)
    erfreq = zeros(5)
    for i in range(len(em2)):
        ncorrect = em2[i,i]
        em2[i,i] = 0
        nwrong = sum(em2[i,:])
        em2[i,:]/=nwrong
        erfreq[i] = nwrong/ncorrect

    

    print em2
    print erfreq

    err = {'error_matrix':error_matrix,
           'error_transitions':em2,
           'error_freqs':erfreq,
           'rights':all_rights,
           'wrongs':all_wrongs,
           'answer':all_cons}
    return err

def view_errors(r):
    print 'MEAN ERROR RATE:'
    print float(np.sum(r.error['wrongs']))/ np.sum(r.error['rights'])

    print 'COMPLETELY CORRECT', float(np.sum(equal(r.error['wrongs'],0))), float(np.sum(equal(r.error['wrongs'],0)))/ len(r.error['wrongs'])*100,'%'
    print 'ONE WRONG', float(np.sum(equal(r.error['wrongs'],1))), 
    print 'TWO WRONG', float(np.sum(equal(r.error['wrongs'],2)))
    print 'MORE THAN TWO', float(np.sum(greater(r.error['wrongs'],2)))

    f = plt.figure(0)
    f.clear()
    ax = f.add_axes([.3,0,1,.7])
    wrongs = r.error['wrongs']
    rights = r.error['rights']
    ax.plot(array(wrongs)/(rights - wrongs)[np.argsort(array(wrongs)/(rights - wrongs))])
    
    ax1 = f.add_axes([0,0,.3,1])
    plt.axis('off')
        
    trans_arrows(r,ax1)



def trans_arrows(r,ax):
    e = r.error['error_transitions']
    
    idx = r.read_idx_dict
    f = plt.figure(1)
    f.clear()
    ni = len(idx.keys())

    ax.set_xlim( [-.5,1.5]) 
    ax.set_ylim([-1,ni] )
    
    #ax.set_aspect('auto')
    #ax.set_autoscale_on(False)
    #ax.set_xlim = [-.5,.5]

    
    for k,v in idx.iteritems():
        probs = e[v,:]
        probs -= min(probs)
        probs /= max(probs)*3

        max_prob = argmax(probs)
        if v == 4: break
        print probs
        for kf,vf in idx.iteritems():
            prob = probs[vf]
            if vf == max_prob:
                alpha = 1
                width = 1
            else:
                alpha = prob
                width = .5

            print prob, v, vf
            p =  patches.ConnectionPatch([0,v],[1,vf],'data','data',
                                         arrowstyle = 'fancy,head_length=0.4,head_width='+str(width)+',tail_width='+str(width),
                                         shrinkA = 50,
                                         shrinkB = 50,
                                   
                                         edgecolor = 'black',
                                         facecolor = 'black',
                                         alpha = alpha)
            ax.add_patch(p)
    for k,v in idx.iteritems():
        ax.annotate(k,[0,v],size = 'xx-large')
        ax.annotate(k,[1,v],size = 'xx-large')
    ax.annotate('likely errors',[.05,.95],xycoords = 'figure fraction',size = 'x-large')
