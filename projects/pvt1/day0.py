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
import Bio.Align as ba
import Bio.Seq as bs
import Bio.SeqRecord as br
import Bio.SeqIO as sio

from numpy import *
import subprocess as spc, re
import StringIO as strio


nt_dict = {
    'A':0,
    'T':1,
    'G':2,
    'C':3,
    'N':4,
    '-':5
    }
nt_rdict = dict([(v,k) for k,v in nt_dict.iteritems()])

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

def run0(spec_ct = 8, **kwargs):

    def setLocusResults(**kwargs):
        spec_ct = kwargs.get('spec_ct')
        bases = (128693265,129266680)
        a0 = fetch_num_ali()
        names = fetch_alinames()
    
        ref = a0[0]
        ali_counts = sum(less(a0,4) * equal(a0, a0[0,:]) ,1)
        names_all = [names[i]  for i in argsort(ali_counts)[::-1]]
        names = names_all[:spec_ct]
        a0 = a0[argsort(ali_counts)[::-1]][:spec_ct]
        ali_counts = sorted(ali_counts)[::-1][:spec_ct]
        
        wl = 150
        n_runs = 500
        locii = {}
        results = {}
        for n_specs in [3, 8]:
             locii[n_specs], results[n_specs] = run_windows(a0,ref, n_specs = n_specs,
                                                            n_runs = n_runs,
                                                            win_len = wl, win_ofs = wl/2,
                                                            spec_names = names)

        return locii, results
    
    return mem.getOrSet(setLocusResults,
                        **mem.rc(kwargs, on_fail = 'compute',
                                 spec_ct = spec_ct))


def check_results(locii, results, n_runs = 400):
    a0 = fetch_num_ali()
    names = fetch_alinames()
            
    f = myplots.fignum(3,(8,8))
    ax = f.add_subplot(211)
    vec = zeros(len(a0[0]))

    xys = {}
    for k in results.keys():
        xys[k] =  array([[l,v['Mean z-score']] for v,l in zip(results[k],locii[k]) if len(v) >= 19],float).T

    raise Exception()
    ax2 = f.add_subplot(111)
    
    ax2.scatter(xys[3][0],xys[3][1])
    #ax2.scatter(xys[8][0],xys[8][1], color = 'red')


    
    f.savefig(myplots.figpath('run0_zscores_{0}runs'.format(n_runs)))

    
    


def run_windows(ali,ref, n_specs = 3,  
                win_len = 150, win_ofs = 75,
                baserng = None,
                spec_names = None,
                window_selection = 'exhaustive',
                n_runs = 2):
    la = shape(ali)[1]
    if baserng == None: baserng = (0,la)
    wvals = []
    ofs = []
    p0 = baserng[0]


    #CHOOSE A LIST OF WINDOW OFFSETS TO INVESTIGATE.
    if window_selection = 'exhaustive':
        while win_len + p0 < baserng[-1]:
            ofs.append(p0)
            p0 += win_ofs
        ofs_choice = ofs
    elif window_selection == 'conserved':
        sum_arr = less(ali,4) * equal(ali,ref)
        ug_ct = less(ali, 4)
        while win_len + p0<baserng:
            wvals.append( sum( sum_arr[:,p0:p0+win_ofs], 1)/( sum(ug_ct[:,p0:p0+win_ofs],1) + .00001))
            p0 = p0 + win_ofs
            ofs.append(p0)
        wvals = array(wvals).T
        tots = sum(wvals[:n_specs],0)
        ofs_choice = argsort(tots)[::-1][:n_runs]
    else:
        raise Exception('window selection type {0} unknown', window_selection)
    
    #INVESTIGATE THE CHOSEN WINDOW SUBSET
    hits = []
    all_vals = []
    for oc in ofs_choice:
        p0 = ofs[oc]
        sub_ali = ali[:n_specs, p0:p0+win_len]
        ali_lets = [[nt_rdict[elt ] for elt in seq] for seq in sub_ali]
        bali = ba.MultipleSeqAlignment([br.SeqRecord(bs.Seq(''.join(let))) for let in ali_lets])
        if spec_names != None:
            for i in range(len(sub_ali)): bali[i].id = spec_names[i]
        child = spc.Popen('RNAz', 
                          stdin = spc.PIPE,
                          stdout = spc.PIPE,
                          stderr = spc.PIPE,
                          shell = True)
        
        aio.write(bali,child.stdin, "clustal")
        child.stdin.close()
        out = child.stdout.read()
        
        pattern = re.compile('^>(?P<name>\S+).*\n'+
                             '(?P<seq>\S+).*\n'+
                             '(?P<struct>\S+)\s*'+
                             '\(\s*(?P<energy>)[^\)]*\)', re.M)
        
        v = list(re.finditer(pattern,out))
        vals = {'seq_props':dict([ (v.groupdict()['name'], v.groupdict()) for v in re.finditer(pattern,out)])}
        vals.update(dict([  [e.strip() for e in elt.split(':')] for elt in out.splitlines() if len(elt.split(':')) ==2]))
        all_vals.append(vals)

    
    return ofs_choice,  all_vals


