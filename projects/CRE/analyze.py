import compbio.config as cfg
from numpy import *
import numpy as np
import cb.utils.memo as mem
import itertools as it
import cb.utils.plots as myplots
import cb.utils.seismic as seismic
import matplotlib.pyplot as plt
from scipy.stats import linregress
import scipy.signal as ss
import cb.utils.sigsmooth as sgs
import subprocess as spc
import cb.utils.colors as mycolors

figtemplate = cfg.dataPath('figs/CRE/{0}.pdf')
promoter_type = 'CRE'

def mut_counts(cons, seqs, name = 'CRE'):
    
    l = len(cons)
    figtitle = '{0}_mut_counts'.format(name)
    f = myplots.fignum(1, (8,8))
    
    ax = f.add_subplot(111)
    counts = zeros(( 4, l))
    lets = ['A','T','G', 'C']
    for pos in range(l):
        counts[:,pos] = [0 if let == cons[pos] else list(seqs[:,pos]).count(let) 
                         for let in lets ]
    seismic.seismic(counts, ax = ax)
    
    f.savefig(figtemplate.format(figtitle))


def position_activities(cons, seqs, activities,
                        show_wt = False):

    
    induction = array([a[0] / a[1] for a in activities])
    l = len(cons)
    
    #wt = [ mean([
    #            val for val in induction[:,i]] if )[nonzero[ for i in range(l)]
    #[for i, let in enumerate(seq)] for seq in seqs.T

    #wt = [[induction[j] 
    #            for j, let in enumerate(seq) if let == cons[i] ] 
    #      for i, seq in enumerate(seqs.T)]
    mut = [[induction[j] 
                 for j, let in enumerate(seq) if let != cons[i] ]
           for i, seq in enumerate(seqs.T)]
    
    percentiles = [10**x for x in range(-5,2) ]+ [40] 
    percentiles = percentiles + [100 -p for p in percentiles]

    f = myplots.fignum(2, (8,8))
    f.clear()
    ax = f.add_subplot(111)
    

    mtiles = [] 
    for i,mt in enumerate( mut):
        if not mt:
            mtiles.append([nan] * len(percentiles))
        else:
            ax.scatter( zeros(len(mt)) + i, \
                            log(mt) + random.random(len(mt))*.1,\
                        2 , color = 'black',alpha = .05)
            mtiles.append([percentile(mt,p) for p in percentiles])

    #ax.plot(arange(len(wmeans)), wmeans, color = 'blue', linewidth = 6,alpha = .2)
    
    mtiles = array(mtiles)
    for mmeans in mtiles.T:
        ax.plot(arange(len(mmeans)), log(mmeans), color = 'red', linewidth = 3, alpha = .6)
    #ax.annotate('$R^2 = {0}$'.format(rsquared),
    #            [1,1], xycoords = 'axes fraction', 
    #            ha = 'right', va = 'top')
    ax.set_xlabel('position mutated (~10 per sequence)')
    ax.set_ylabel('log induction (induced expr/uninduced)')
    ax.set_title('Induction ratios for sequences mutated at points')
    
    figtitle = 'single_position_percentiles'
    f.savefig(figtemplate.format(figtitle))

def position_interactions(cons, seqs, activities):
    'POSITION INTERACTIONS PA/PB on X Y axes.'
    l = len(cons)

    
    mut = [[induction[j] 
                 for j, let in enumerate(seq) if let != cons[i] ]
           for i, seq in enumerate(seqs.T)]
    mut_pairs = ''
    
    raise Exception()


def meme_finder():
    'Run meme on high induction random mutations of the sequence.'

    pass

def transfac_finder():
    'Then maybe run the thing through transfac. The transcription factor database.'

    pass




#ENERGY ACTIVITY STUFF.
def site_energy_deltas(showtype = 'first_part_energies',
                       muts_allowed = None,
                       mkey = None,
                       figtitle = 'sed',
                       smoothing = None,
                       sub_means = True,
                       en_type = 'double',
                       induction_type = 'ratio'):
    cre, cre_rndvals, keys = get_mutants()
    mut_inds = site_mut_inds()
    
    mean_induction = get_mean_induction()
    xs, ys, rs, cs , pas, mas, was, bads = [[] for i in range(8)]

    cons = get_cons()
    l = len(cons)
    energies = get_energies()
    enfun = lambda x : sum([energies[''.join(x[pos:pos+2])][2] for pos in range(len(x) -1)])
    
    flip = 0
    
    figtitle = figtitle + '_' + showtype
    if not mkey == None:
        muts_allowed = set(get_motif_dicts()[mkey])
        figtitle = figtitle + '_{0}'.format(mkey)
    if not smoothing  == None:
        figtitle = figtitle + '_sm={0}'.format(smoothing)
        
    if showtype == 'first_part_energies':
        if en_type == 'double':
            ens = get_energies()
            gibbs = dict([(k,v[2]) for k,v in ens.iteritems()])
        else:
            gibbs = get_sing_energies()
        figtitle = figtitle + '_etype={0}'.format(en_type)
        enlen = len(gibbs.keys()[0])
        doubles = gibbs.keys()
        keysrt = dict([(doubles[i], elt) for i,elt in enumerate( argsort(argsort(gibbs.values())))])
        dub_muts = [[[] for j in range(len(keysrt)) ] for i in range(len(keysrt))]

    figtitle = figtitle + '_ind=' + induction_type

    for site, muts in enumerate(mut_inds):
        if muts_allowed != None: muts = list(muts_allowed.intersection(muts))
        mut_avg = mean( cre_rndvals[muts,0] / cre_rndvals[muts,1])
        
        trip_rng = position_triplet(site)
        wt_seq = array([cons[s] for s in trip_rng])

        mut_seqs = cre[muts][:,trip_rng]
        
        if induction_type == 'ratio':  mut_inductions =(cre_rndvals[muts,0] / cre_rndvals[muts,1])
        elif induction_type == 'on': mut_inductions =  cre_rndvals[muts,0]
        elif induction_type == 'off': mut_inductions = cre_rndvals[muts,1]

        wt_e = enfun(wt_seq)
        mut_es = [enfun(seq) for seq in mut_seqs]
        
        if showtype == 'show_avgs':
            
            plus_avg = mean(mut_inductions[nonzero(greater(mut_es, wt_e))[0]])
            minus_avg = mean(mut_inductions[nonzero(less_equal(mut_es, wt_e))[0]])
            cs.extend(['red', 'blue'])
            xs.extend([site] * 2)
            ys.extend([plus_avg, minus_avg])
            rs.extend([[30] * 2])
            
        elif showtype == 'show_all':
                                                    
            cs.extend(['red' if e > wt_e else 'blue' for e in mut_es])
            xs.extend([site] * len(mut_es))
            ys.extend(mut_inductions + np.random.rand(len(mut_es)) * .1)
            rs.extend([10] * len(mut_es))
    
        elif showtype == 'smoothed_avg':
            
            plus_avg = mean(mut_inductions[nonzero(greater(mut_es, wt_e))[0]])
            minus_avg = mean(mut_inductions[nonzero(less_equal(mut_es, wt_e))[0]])
            whole_avg = mean(mut_inductions)

            if isnan(plus_avg): 
                bads.append(site)
                plus_avg = 1
            if isnan(minus_avg):
                bads.append(site)
                minus_avg = 1
            if isnan(whole_avg): 
                bads.append(site)
                whole_avg = 1

            pas.append(plus_avg)
            mas.append(minus_avg)
            was.append(whole_avg)


            cs.extend(['red', 'blue'])
            xs.extend([site] * 2)
            ys.extend([plus_avg, minus_avg])
            rs.extend([[30] * 2])
                        
        elif showtype == 'first_part_energies':
            for idx, trip in enumerate(mut_seqs):
                if len(trip) == 2 : continue
                if site < 1: continue
                if site > 28: continue

                for pos in range(len(trip) - enlen):
                    midx = keysrt[''.join(trip[pos:pos+enlen])]
                    cidx = keysrt[''.join(wt_seq[pos:pos+enlen])]
                    dub_muts[cidx][midx].append( mut_inductions[idx] )
            
    
    f = myplots.fignum(3,(8,8))
    ax = f.add_subplot(111)

    

    if showtype in ['show_avgs', 'show_all']: 
        ydat = np.log(ys)
        xdat = array(xs)
        cs = array(cs)
        rs = array(rs)
        ax.scatter(xdat, ydat, array(rs), c= cs, alpha = .2)
        
        xlim = [min(xs), max(xs)]
        ylim = [min(ydat), max(ydat)]
        ax.set_title('Average induction of site mutants (red mutants have stronger pairing)')
        ax.set_ylabel('fold induction')
        ax.set_xlabel('mutation site')
            

    elif showtype == 'smoothed_avg':
        pas =log( array(pas))
        mas =log( array(mas))
        was =log(array(was))
        if sub_means: subvec = was
        else: subvec = zeros(len(was))
        if smoothing == None:
            plus = pas - subvec
            minus =mas - subvec 
        else:
            plus =sgs.smooth(pas  -subvec, smoothing)
            minus =sgs.smooth(mas -subvec, smoothing)
    
        goods = ones(len(plus))
        goods[[b for b in bads]] = 0
        plus = plus * goods
        minus = minus * goods

        ydat = plus
        xdat = arange(len(plus))

        ax.plot(plus , color = 'white', linewidth = 5)
        ax.plot(minus, color = 'white', linewidth = 5)
        ax.plot(plus, color = 'red', linewidth = 3)
        ax.plot(minus, color = 'blue', linewidth = 3)
        ax.fill_between(arange(l), plus, minus, where = greater_equal(minus, plus),
                        color = 'blue', alpha = .5, interpolate = True)
        ax.fill_between(arange(l), minus, plus, where = less_equal(minus, plus),
                        color = 'red', alpha = .5, interpolate = True)
            
        #ax.set_xlim([min(xdat),max(xdat)]); ax.set_ylim([min(ydat), max(ydat)])
        ax.set_title('Average induction of site mutants (red mutants have stronger pairing)')
        ax.set_ylabel('fold induction')
        ax.set_xlabel('mutation site')
    elif showtype == 'first_part_energies':
        f = myplots.fignum(3,(8,4))
        f.clear()
        ax = f.add_subplot(121)
        img =array( [[mean(log(inductions)) if sum(inductions) != 0 else nan for inductions in row] for  row in dub_muts])
        img[isinf(img)] = min(img[isfinite(img)]) 
        ax.imshow(img,
                  interpolation = 'nearest',
                  cmap = plt.get_cmap('OrRd'))
        ax2 = f.add_subplot(122)
        ax2.imshow([[len(log(inductions)) for inductions in row] for row in dub_muts],
                  interpolation = 'nearest',
                  cmap = plt.get_cmap('OrRd'))
        for a in [ax,ax2]:
            a.set_xticks(keysrt.values())
            a.set_xticklabels(keysrt.keys())
            a.set_yticks(keysrt.values())
            a.set_yticklabels(keysrt.keys())
        ax.set_title('mean induction change')
        ax2.set_title('transition counts')
        ax.set_ylabel('wildtype base')
        ax.set_xlabel('mut base')

    else:
        raise Exception()

    if showtype in ['show_avgs', 'show_all', 'smoothed_avg']:
        l0 = [ax.get_ylim(), ax.get_xlim()]
        
        for rng in cre_rngs():
            ax.plot(rng,[0]*2, linewidth = 6, color = 'black') 
            ax.fill_betweenx(ax.get_ylim(), [rng[0]] * 2, [rng[1]] * 2, alpha = .2, color = 'black')
            for r in rng:
                ax.vlines(r,*ax.get_ylim(),alpha = .75)
        ax.set_ylim(l0[0]); ax.set_xlim(l0[1])
        


    ax.annotate(figtitle, [0,0],
                xycoords ='figure fraction',
                va = 'bottom', ha = 'left')

    f.savefig(figtemplate.format(figtitle))


#UTILS
def load_motifs():
    mfile = open(cfg.dataPath('motifs/vert_tfs.txt'))
    mdicts = {}
    cur_lets = None
    for l in mfile.xreadlines():
        if l[0] == '>':
            if cur_lets != None:
                mdicts[motif_key]['lets'] = cur_lets
                mdicts[motif_key]['pwm'] = array(pwm)
            motif_key = l.split(' ')[0][1:]
            mdicts[motif_key] = {}
            cur_lets = []
            pwm = []
        else:
            cur_let, cur_pwm = l[0:1], l[2:].split(' ')
            cur_lets.append(cur_let)
            pwm.append(cur_pwm)

    raise Exception()
            
def write_seqs_to_motifs():
    cre, rnd, keys = get_mutants()
    cons = get_cons()
    
    contents = ''
    for i, c in enumerate(cre):
        k = keys[i]
        name = k
        contents +=  '\n'.join(['A {0} 1 {1}'.format(k,len(cons)),
                           '>{0}'.format(promoter_type),
                           ''.join(c.lower()),'\n'])
    outfile = open(cfg.dataPath('CRE/{0}_for_motifs.txt'.format(promoter_type)),'w')
    outfile.write(contents)

def get_motifs(**kwargs):
    def set_motifs(**kwargs):
        fpath = cfg.dataPath('CRE/{0}_for_motifs.txt'.format(promoter_type))
        cmd = 'motif-match -n 1 -m {0} -V 1'.format(fpath)
        out = spc.Popen(cmd, shell = True, stdout = spc.PIPE)
        comm = out.communicate()
        return comm
        
    return mem.getOrSet(set_motifs, **mem.rc(kwargs,
                                             on_fail = 'compute',
                                             register = promoter_type))

def get_mutants(**kwargs):
    if promoter_type == 'CRE':
        return getCRE(**kwargs)
    else: 
        raise Exception()
        

def nt_ids():
    return {'A': 0,
            'T': 1,
            'G': 2,
            'C': 3}
def id_nts():
   return array(['A','T','G','C'])

def get_energies():

    return  {'AA':[8.0,21.9,1.2],
           'TT':[8.0,21.9,1.2],
           'AT':[5.6,15.2, .9],
           'TA':[6.6,18.4, .9],
           'CA':[8.2,21.0,1.7],
           'TG':[8.2,21.0,1.7],
           'CT':[6.6,16.4,1.5],
           'AG':[6.6,16.4,1.5],
           'GA':[8.8,23.5,1.5],
           'TC':[8.8,23.5,1.5],
           'GT':[9.4,25.5,1.5],
           'AC':[9.4,25.5,1.5],
           'CG':[11.8,29.0,2.8],
           'GC':[10.5,26.4,2.3],
           'GG':[10.9,28.4,2.1],
           'CC':[10.9,28.4,2.1]}
    
def get_sing_energies():
    e0 = get_energies()
    eout = dict([(k , mean([ elt[2] for dub,elt in e0.iteritems() if k in dub])) for k in ['A','T','G','C']])
    return  eout
    


def cre_rngs():
    return [(11,19),
            (37,42),
            (47,52),
            (69,77)]
def cre_masks(pad = 2):
    cons = get_cons()
    l = len(cons)
    masks = zeros((4,l))
    for i,c in enumerate( cre_rngs()):
        c=array(c) + [-pad , pad]
        masks[i][arange(*c) ] = 1
    return masks

def get_motif_dicts(pad = 2, **kwargs):
    def set_motif_dicts(**kwargs):
        masks = cre_masks(kwargs.get('pad'))
        out = {}
        cons = [nt_ids()[let] for let in get_cons()]
        for j, seq in enumerate(get_num_seqs()):
            key = tuple([ i for i , mask in enumerate( masks ) 
                          if sum(not_equal(seq,cons) * mask) != 0 ]) 
            if not out.has_key(key): out[key] = []
            out[key].append(j)
        return out
    return mem.getOrSet(set_motif_dicts,
                        **mem.rc(kwargs,
                                 pad = pad,
                                 register = '{0}_{1}'.format(promoter_type, pad),
                                 on_fail = 'compute'))

def filters(name, num = 3, nums = (1,2)):
    seqs = get_num_seqs()
    cons = array([nt_ids()[let] for let in get_cons()])
    l = len(cons)
    if name == 'no_motif':
        '''No motif mutations'''
        mask = zeros(l)
        for r in cre_rngs(): mask[arange(*r)] = 1
        out =  equal(0, sum(not_equal(seqs, cons) * mask,1))
    
    elif name == 'with_motifnum':
        '''Demand at least on mutation in the specified motif'''
        mask = cre_masks()[num]
        out = greater(  sum(not_equal(seqs, cons) * mask,1), 0)
        
    elif name == 'without_motifnum':
        '''No mutation in the specified motif'''
        mask = cre_masks()[num]
        out = equal( sum(not_equal(seqs, cons) * mask, 1),0)
        
    elif name == 'with_motifnums':
        '''Demand at least one mutation in each of the specified motif'''
        mask = sum([masks[n] for n in nums], 0) 
        out = greater(  sum(not_equal(seqs, cons) * mask,1), 0)

    elif name == 'only_motifnum':
        masks = cre_masks()
        mask1 = greater( sum(not_equal(seqs, cons) * masks[num],1), 0)
        mask2 = equal(  sum(not_equal(seqs, cons) * \
                                sum([masks[i] 
                                     for i in range(len(masks))
                                     if not i == num ]
                                    ,0), 1),0)
        out = mask1 * mask2 
    elif name == 'only_motifnums':
        masks = cre_masks()

        mask1 = not_equal(  sum(not_equal(seqs, cons) * \
                                    sum([masks[i] 
                                     for i in range(len(masks))
                                     if i in nums ]
                                    ,0), 1),0)

        mask2 = equal(  sum(not_equal(seqs, cons) * \
                                sum([masks[i] 
                                     for i in range(len(masks))
                                     if i not in nums ]
                                    ,0), 1),0)
        out = mask1 * mask2
        
    
    return out
def get_num_seqs(**kwargs):
    def set_num_seqs(**kwargs):
        ntdict = nt_ids()
        cre, cre_rndvals, keys = get_mutants()
        return array([[ntdict[let] for let in seq] for seq in cre])
    return mem.getOrSet(set_num_seqs, 
                        **mem.rc(kwargs,
                                 register = promoter_type,
                                 on_fail = 'compute'))
                     

def position_triplet(idx):
    cons = get_cons()
    l = len(cons)
    if idx == 0: trip_rng = arange(0,2)
    elif idx == l-1: trip_rng = arange(l-2,l)
    else: trip_rng = arange(idx -2, idx+1)    
    return trip_rng

def get_trip_muts(idx):
    cre, cre_rndvals, keys = get_mutants()
    seqs = get_num_seqs()
    cons_let = get_cons()
    ntdict = nt_ids()
    cons = array([ntdict[let] for let in cons_let])
    l = len(cons)
    trip_rng = position_triplet(idx)
    
    muts = nonzero( not_equal(np.sum( seqs[:,trip_rng ] - cons[trip_rng], 1),0))[0]
    return muts
    

def site_mut_inds(**kwargs):
    def set_site_muts(**kwargs):
        l = len(get_cons())
        site_muts = [set( get_trip_muts(idx) ) for idx in range(l)]
        return site_muts
    return mem.getOrSet(set_site_muts, **mem.rc(kwargs,
                                                register = promoter_type,
                                                on_fail = 'compute'))
def get_mean_induction(**kwargs):
    def set_mind(**kwargs):
        cre, cre_rndvals, keys = get_mutants()
        return mean(cre_rndvals[:,0])/ mean(cre_rndvals[:,1])
    return mem.getOrSet(set_mind, **mem.rc(kwargs, 
                                           register = promoter_type,
                                           on_fail = 'compute'))

def getCRE(**kwargs):
    def setCRE(**kwargs):        
        cre_des = open(cfg.dataPath('CRE/27k/CRE_Randomization_Design.txt'))
        cre_rnd = open(cfg.dataPath('CRE/27k/CRE_Randomization.dat'))
        cre_rnd = cre_rnd.readlines()
        cre_des = cre_des.readlines()

        cre_rndvals = [[elt.strip() for elt in line.split('\t')] for line in cre_rnd[1:]]
        cre_seqs = [[elt.strip() for elt in line.split('\t')] for line in cre_des]

        cre_rndvals = dict([(e[0], e[1:]) for e in cre_rndvals])
        cre_seqs = dict([[e[0],e[1]] for e in cre_seqs])

        keys = list(set(cre_rndvals.keys()).intersection(cre_seqs.keys()))
        
        cre = array([list(cre_seqs[k]) for k in keys])
        cre_rndvals = array([array(cre_rndvals[k], float) for k in keys])
        return cre, cre_rndvals, keys
    return mem.getOrSet(setCRE,
                        **mem.rc(kwargs,
                                 register = promoter_type,
                                 on_fail = 'compute'))


def get_cons(**kwargs):
    def consensus_seq(seqs):
        return [ sorted([(k,list(g)) for k, g in it.groupby(sorted(c)) ], key = lambda x: len(x[1]))[-1][0] for c in seqs.T]
    def set_cons(**kwargs):
        cre, cre_rndvals, keys = get_mutants(**mem.sr(kwargs))
        cons = consensus_seq(cre[::100])
        return cons
    cons = mem.getOrSet(set_cons, **mem.rc(kwargs,
                                           register = promoter_type,
                                           on_fail = 'compute'))
    return cons


