#Get a list of SNPs for all of the hapmap populations.
import re, os, gzip, itertools as it
import compbio.config as cfg
from numpy import *
import scipy.linalg as la

import matplotlib.pyplot as plt
import cb.utils.plots as myplots
import cb.utils.colors as mycolors

import cb.utils.elpsfit as elpsfit
from cb.utils.seismic import seismic
'''
Workflow:

Calling individuals = hapmap.individuals() returns a dict of genotypes keyed
by hapmap sample ID, each entry having a list of SNP genotypes in alphabetic 
form.

To go from an alphabet to an array, call hapmap.indy_array(individuals) to
return an integer array having values set according to the base map (ATGCN).

Additional operations include computing principal components of variation
and projecting onto principal components. 

Clustering to find locii easily describable in terms of a small number of
clustered haplotypes for the purpose of selectivity computation is not yet
implemented. 

Should figure out how the linkage disequilibrium is traditionally
computed in order to skip having to compute this myself. Also, note the
definition of F_st to hightlight regions having anomolously distributed SNP
haplotypes. Not implemented.

'''

def region_abbrevs(collection = 'hapmap'):
    if collection == 'hapmap':
        region_text = '''     Maasai in Kinyawa, Kenya (MKK)
	Yoruba in Ibadan, Nigeria (YRI)
	Han Chinese in Beijing, China (CHB)
	Japanese in Tokyo, Japan (JPT)
	Luhya in Webuye, Kenya (LWK)
	Chinese in Metropolitan Denver, CO, USA (CHD)
	Gujarati Indians in Houston, TX, USA (GIH)
	Toscani in Italia (TSI)
	Mexican Ancestry in Los Angeles, CA, USA (MXL)
	African Ancestry in SW USA (ASW)'''

    region_list =[re.compile( '(?P<description>.*)\((?P<abbrev>.*)\)' )\
                      .search(l).groupdict()
                  for l in region_text.splitlines()]
    return dict([(e['abbrev'],e['description'].strip()) for e in region_list])


def polyfile(chr = 10, grp = 'ASW', **kwargs):
    root = cfg.dataPath('hapmap/phase3/polymorphic')
    fname = os.path.join(root,
                         'genotypes_chr{0}_{1}_phase3.2_nr.b36_fwd.txt.gz'\
                             .format(chr,grp));
    contents = gzip.open(fname).readlines()                    
    return contents

def consfiles(chr = 10, maxlines = 1000, output = 'list'):
    root = cfg.dataPath('hapmap/phase3/consensus')
    chr_regex = re.compile('_chr{0}_'.format(chr))
    all_content = []
    fnames = []
    for r,d,files in os.walk(root):
        for f in files:
            if chr_regex.search(f):
                path = os.path.join(r,f)
                fopen =gzip.open(path)
                all_content.append(fopen.readlines()[:maxlines])
                fopen.close()
                fnames.append(path)
    if output == 'list':
        return all_content
    else:
        return dict([(fnames[i], c) for i, c in enumerate(all_content)])

def polyfiles(chr = 10, maxlines = 1000, output = 'list'):
    root = cfg.dataPath('hapmap/phase3/polymorphic')
    chr_regex = re.compile('_chr{0}_'.format(chr))
    all_content = []
    fnames = []
    for r,d,files in os.walk(root):
        for f in files:
            if chr_regex.search(f):
                path = os.path.join(r,f)
                fopen =gzip.open(path)
                all_content.append(fopen.readlines()[:maxlines])
                fopen.close()
                fnames.append(path)
    if output == 'list':
        return all_content
    else:
        return dict([(fnames[i], c) for i, c in enumerate(all_content)])

def indy_info():
    infos = {}
    all_fcontents = polyfiles(maxlines = 1, output = 'dict')
    for fname, fcontents in all_fcontents.iteritems():
       cols = fcontents[0].split(' ')
       nattr = 11
       attr_cols = cols[:nattr]
       indy_cols = [c.strip() for c in cols[nattr:]]
       region_symb = re.compile('chr\d+_([A-Z]*)_').search(fname).groups()[0]
       region_desc = region_abbrevs().get(region_symb, 'description unavailable')
       for k in indy_cols:
           infos[k] = {'fname': fname,
                       'region_desc':region_desc,
                       'region_symb':region_symb}

    keys = infos.keys()
    srt = argsort(keys)
    for i,idx in enumerate(srt):
        infos[keys[idx]]['sortorder'] = i
    
    return infos
                       




       

def individual_snp_dicts(**kwargs):
    all_fcontents =polyfiles(**kwargs)
    all_individuals = []
    all_snps = []
    
    for fcontents in all_fcontents:
       cols = fcontents[0].split(' ')
       nattr = 11
       attr_cols = cols[:nattr]
       indy_cols = [c.strip() for c in cols[nattr:]]
       
       individuals = {}
       snps = []
       
       for i in indy_cols:
           individuals[i] = {}
           
       for l in fcontents[1:]:
           snp =dict(zip(attr_cols, l.split(' ')[:nattr]))
           snps.append(snp)
           rs = snp['rs#']

           ls =  l.split(' ')
           for i,k in enumerate(indy_cols):
               individuals[k][rs] =ls[i + nattr]
       
       all_individuals.append(individuals)
       all_snps.append(snps)
       
       
    individuals = {}
    for item in all_individuals:
        individuals.update(item)
    snps = {}
    for snplist in all_snps:
        snps.update(dict([(elt['rs#'],elt ) for elt in snplist]))

    return individuals, snps

      
       

def full_snp_dicts(**kwargs):
    cons,csnps = consensus_snp_dicts()
    indy,isnps = individual_snp_dicts()
    out = {}
    for k in indy.keys():
        out[k] = {}
        out[k].update(cons[k])
        out[k].update(indy[k])

    snps = {}
    snps.update(csnps)
    snps.update(isnps)
    return out, snps


def consensus_snp_dicts(**kwargs):
    all_fcontents =consfiles(**kwargs)
    all_individuals = []
    all_snps = []
    
    for fcontents in all_fcontents:
       cols = fcontents[0].split(' ')
       nattr = 11
       attr_cols = cols[:nattr]
       indy_cols = [c.strip() for c in cols[nattr:]]
       
       individuals = {}
       snps = []
       
       for i in indy_cols:
           individuals[i] = {}
           
       for l in fcontents[1:]:
           snp =dict(zip(attr_cols, l.split(' ')[:nattr]))
           snps.append(snp)
           rs = snp['rs#']

           ls =  l.split(' ')
           for i,k in enumerate(indy_cols):
               individuals[k][rs] =ls[i + nattr]
       
       all_individuals.append(individuals)
       all_snps.append(snps)
       
       
    individuals = {}
    for item in all_individuals:
        individuals.update(item)
    snps = {}
    for snplist in all_snps:
        snps.update(dict([(elt['rs#'],elt ) for elt in snplist]))

    return individuals, snps



def indy_array(individuals , linear = True, strands = False):
    '''
    Get a olymorphism matrix describing each of the individuals.

kwargs:
linear:  flatten the last dimension of the array return vectors for each individual.
         default is true, if false, retain the [ni x l x 5] shape giving matrices per
         individual.

returns:
           individual SNP genotypes according to the base map provided.
           [ni x (l*5)] or [ni x l x 5] depending on the value of the linear kwd
'''
    snp_keys = set(list(it.chain(*[d.keys() for d in individuals.values()])))
    

    base_map = {'A':0,
                'T':1,
                'G':2,
                'C':3,
                'N':4}

    ikeys = sorted(individuals.keys())
    snp_ids = sorted(snp_keys)
    l = len(snp_ids)

    if strands:
        out = zeros((len(ikeys)*2,l,5))
    else:
        out = zeros((len(ikeys),l,5))


    for i,k in enumerate(sorted(ikeys)):
        for j, rs in enumerate(snp_ids):
            if strands:
                out[2*i,j,base_map[individuals[k].get(rs,'NN')[0]]] = 1 
                out[2*i + 1,j, base_map[individuals[k].get(rs,'NN')[1]]]= 1
            else:
                out[i,j,base_map[individuals[k].get(rs,'NN')[0]]] = 1 

    if not linear: return out
    s = shape(out)
    l = s[1]
    ni = s[0]
    nb = s[2]
    vecs= reshape(out,(ni, l * nb))
    return vecs

def cc(indy_arr):
    '''
    output of indy_arr must be set to linear.
'''
    return corrcoef(indy_arr)


def cc_svd(cc, indy_arr):
    U,s,Vh = la.svd(cc)
    pcs = U.T

    gene_pcs = dot(indy_arr.T, pcs[:5].T)
    return pcs, gene_pcs, s


def indy_regions(indy_arr, indy_info):
    regions_sorted = sorted(set([x['region_symb'] for x in indy_info.values()]))
    regions_map = dict([(k, i) for i , k in enumerate(regions_sorted)])

    nr = len(regions_sorted)
    ni = len(indy_info.keys())

    indy_regions = zeros(ni)
    for k,v in indy_info.iteritems():
        row = v['sortorder']
        region = v['region_symb']
        indy_regions[row] = regions_map[region]
        
    return array(indy_regions, int) 

def snp_count_plot(indy_arr, indy_info):
    regions = indy_regions(indy_arr, indy_info)
    f = myplots.fignum(3, (8,8))
    ax = f.add_subplot(111)
    for row in indy_arr.T[4::5][:50]:
        rsub = array(regions[:500],float) / max(regions[:500])
        inds = argsort(rsub)
        ax.plot(row[:500][inds] + random.rand(500) * .1)
        ax.plot(rsub[inds], linewidth = 5)
    return

def snp_counts(indy_arr, indy_info):
    
    regions = indy_regions(indy_arr, indy_info)
    
    ct = mycolors.getct(len(regions))
    skip = 1
    ofs = 4
    
    f = myplots.fignum(4, (8,8))
    ax = f.add_subplot(111)
    
    n_snps = 20
    rset = set(regions)
    rcounts = zeros((len(rset), n_snps))

    xs = []
    ys = []
    cs = []
    rs = []
    
    for i, snp in enumerate(indy_arr.T[4::5][:50]):
        rsub = array(regions[::100],float) / max(regions[:20])
        inds = argsort(rsub)
        ys.extend([i] *len(rsub))
        rs.extend([10 + 30 * (1 - snp[:2:100])])
        xs.extend(rsub)
        
        print i
        
    ax.scatter(xs,ys,rs)


    f.savefig(myplots.figpath('regional_snp_counts_first.pdf'))
    
    return


def map_pca(indy_arr, indy_info, pca = None,
            metagroup = 'all',
            **kwargs):
    '''
    output of indy_arr must be set to linear.
'''
    
    metagroups = region_groups(indy_info)
    if metagroup == 'all':
        inds_allowed = arange(len(indy_arr))
    else:
        inds_allowed = metagroups[metagroup]['indy_idxs']

    #fetch some global variables
    regions = array(indy_regions(indy_arr, indy_info))
    iis = [x[1] for x in sorted([(k,v) for k,v in indy_info.iteritems()],
                           key = lambda x:x[0])]
    #and then extract the metaregion
    indy_arr = indy_arr[inds_allowed]
    regions = regions[inds_allowed]
    iis = [iis[ia] for ia in inds_allowed]
    
    #compute pca if none is given
    if pca == None:
        corr = cc(indy_arr)
        indy_pcs, gene_pcs, s = cc_svd(corr, indy_arr)
        pca = gene_pcs

    #and from pca, coordinates
    xys = dot(indy_arr, pca[:,:2]).T


    #then plot
    ct = array(mycolors.getct(max(regions)+1))
    r_abbrevs = region_abbrevs()

    f = myplots.fignum(1,(14,8))
    ax = f.add_subplot(111)
    ax.set_title('top two principle component projection of individual genotypes')

    for k,g in it.groupby(sorted(enumerate(regions),
                                 key = lambda x: x[1]),
                          key = lambda x: x[1]):
        grp = list(g)
        e0 = grp[0][0]
        i0 = iis[e0]
        ax.scatter(xys[0,[e[0] for e in grp]], 
                   xys[1,[e[0] for e in grp]]
                   , color = ct[k],
                   label = '{0} - {1}'.format(i0['region_symb'],
                                              i0['region_desc']))

    #plot regionwide ellipses

    c_imap = dict([(v,k) for k,v in enumerate(set(regions))])
    c_rmap = dict([(v,k) for k,v in c_imap.iteritems()])
    clusters = array([c_imap[r] for r in regions])
    csort = argsort(clusters)

    elps, infos = elpsfit.cluster_ellipses(xys[:,csort],clusters[csort])
    for i, e in enumerate(elps): 
        e.set_alpha(.75)
        e.set_facecolor(ct[c_rmap[i]])
        e.set_edgecolor('none') 
        maj = e.width
        minor = e.height
        e.height = maj * 2
        e.width = minor * 2
    
        ax.add_patch(e)

    ax.legend()

    f.savefig(myplots.figpath('pca_grouping={0}.pdf'.format(metagroup)))
    return

def region_groups(indy_info):
    grps =  {'africa':{'symbs':['ASW',
                              'LWK',
                              'MKK',
                              'YRI']},
             'other':{'symbs':['CEU',
                             'CHB',
                             'CHD',
                             'GIH',
                             'JPT',
                             'MEX',
                             'TSI']}}
    for k,v in grps.iteritems():
        v['indy_idxs'] = [indy['sortorder']
                          for indy in indy_info.values()
                          if indy['region_symb'] in v['symbs']]
    return grps
    

def cluster_range(indy_arr):
    import rpy2.robjects as robjects
    from rpy2.rlike.container import TaggedList
    from rpy2.robjects.packages import importr
    r = robjects.r
    hcv = r.hclust(r.dist(mb))
