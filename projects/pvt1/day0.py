'''
Load up the part of the 32 mammals genome corresponding to PvT1.

Check coverage.

Save a CLW file.

Run RNAz on windows.

Then run RNAz for long range interactions.

'''
import cb.config as cfg
import cb.utils.plots as myplots
import cb.utils.memo as mem
import Bio.AlignIO as aio
from numpy import *


nt_dict = {
    'A':0,
    'T':1,
    'G':2,
    'C':3,
    'N':4,
    '-':5
    }

def fetch_num_ali():
    def setFetchNumAli(**kwargs):
        fa_file = cfg.dataPath('pvt1/pvt1.fa')

    
        ali = aio.parse(open(fa_file), 'fasta')
        a_ali = ali.next()
    
        a0  =  [[nt_dict[elt] for elt in a.seq.upper()] for a in a_ali]
        a0_num = array(a0, byte)        
        return a0_num
    return mem.getOrSet(setFetchNumAli,
                        on_fail = 'compute')
def fetch_alinames():
    def setFetchAliNames(**kwargs):
        fa_file = cfg.dataPath('pvt1/pvt1.fa')

    
        ali = aio.parse(open(fa_file), 'fasta')
        a0 = ali.next()
        return [a.id for a in a0]

    return mem.getOrSet(setFetchAliNames,
                        on_fail = 'compute')    

def run0(spec_ct = 8):
    bases = (128693265,129266680)
    a0 = fetch_num_ali()
    names = fetch_alinames()

    ref = a0[0]
    ali_counts = sum(less(a0,4) * equal(a0, a0[0,:]) ,1)
    names_all = [names[i]  for i in argsort(ali_counts)[::-1]]
    names = names_all[:spec_ct]
    a0 = a0[argsort(ali_counts)[::-1]][:spec_ct]
    ali_counts = sorted(ali_counts)[::-1][:spec_ct]
    
    wl = 100
    win_counts = run_windows(a0,ref, win_len = wl, win_ofs = wl/2)

    f = myplots.fignum(3,(8,8))
    ax = f.add_subplot(211)
    vec = zeros(len(a0[0]))
    ax.plot(ali_counts)
    
    ax2 = f.add_subplot(212)
    ax2.imshow(win_counts,
               aspect = 'auto',
               interpolation = 'nearest')


    ax.set_xticklabels(names, rotation = 45,
                       ha = 'right', va = 'top')
    ax.set_xticks(range(len(names)))

    ax2.set_yticklabels(names, rotation = 45,
                       ha = 'right', va = 'top')
    ax2.set_yticks(range(len(names)))
    f.savefig(myplots.figpath('run0_alignment_hits'))
    
    



    


def run_windows(ali,ref,  win_len = 100, win_ofs = 50):
    la = shape(ali)[1]
    p0 = 0
    wvals = []
    
    sum_arr = less(ali,4) * equal(ali,ref)
    ug_ct = less(ali, 4)
    while win_len + p0< la:
        wvals.append( sum( sum_arr[:,p0:p0+win_ofs], 1)/( sum(ug_ct[:,p0:p0+win_ofs],1) + .00001))
        p0 = p0 + win_ofs
    return array(wvals).T
                    


