import cb.config as cfg
import numpy as np
from numpy import *
import os, itertools as it
import cb.utils.memo as mem

import Bio.SeqIO as sio
import gzip
import cb.external.GFF as GFF
import re

'''
parse all the data required to produce a 
worm chip network
'''

default_atype = 'wormtile'

def tf_chip_peaks(**kwargs):
    def setTf_Chip_Peaks(**kwargs):
        root = cfg.dataPath('wormchip')
        files = [os.path.join(root, f) for f in os.listdir(root)]
        out = {}
        for f in files:
            fopen= open(f)
            data = [l for l in fopen.readlines() if not l[0] == '#']
            
            out[os.path.basename(f)] = \
                [dict(zip(['chr', 'meth', 'type',
                  'start','end','score',
                  'blank','blank2','qValue' ], 
                         l.strip().split('\t')))
                 for l in data]
            vlens = [len(e) for e in out[os.path.basename(f)]]

        for k,v in out.iteritems():
            for d in v:
                d['start'] = int(d['start'])
                d['end'] = int(d['end'])
                d['score'] = float(d['score'])
                d['qValue'] = float(d['qValue'].split('=')[1])
        return out
    return mem.getOrSet(setTf_Chip_Peaks, 
                        **mem.rc(kwargs,
                                 hardcopy = True))
#deprecated... use get_assay info instead
def tf_chips(**kwargs):
    def set_tf_chips(**kwargs):
       #peaks and keys 
       assay_peaks = tf_chip_peaks()
       pkeys = assay_peaks.keys()
       
       #index and uniquify tfs
       ptfs = dict([(k,k[:k.index(':')]) for k in pkeys])
       tfnames = set(ptfs.values())
       tfkeys= dict([(name, [k for k,v in ptfs.iteritems()
                             if v == name])
                     for name in tfnames])       

       chrmap = dict(zip(['I','II', 'III','IV', 'V', 'X'],
                         chromosome_names()))
       
       all_tss = get_tss()
       all_genes = parse_genes()
       
       bind_infos = {}
       #loop through tfs, assays, and finally peaks
       for tfname in tfnames:
           print tfname
           bind_infos[tfname] = {}
           for ekey in tfkeys[tfname]:
             print 'n_exps = {0}'.format(np.sum([len(v) 
                                              for v in bind_infos.values()]))
             exp = assay_peaks[ekey]
             bind_infos[tfname][ekey] = []
             for e in exp:
               start = e['start']
               end = e['end']
               c = e['chr']

               #if we have non chromosomal DNA, skip it!
               if c in chrmap.keys():
                   crkey = chrmap[c]
               else: continue
               genes = all_genes[crkey]
               tss = all_tss[crkey]
               
               #grab indexes in the list of upstream and ds genes
               fdownstream_idx =searchsorted(tss['fwd_tss'],(start+end) /2)
               rdownstream_idx =searchsorted(tss['rev_tss'],(start+end) /2, 'right')-1
               
               #handle the case where the factor hits past the 
               #last gene
               fup_ofs = -1
               if fdownstream_idx == len(tss['fwd_genes']):
                   fdownstream_idx = fdownstream_idx -1
                   fup_ofs = 0
               
               rup_ofs = 1
               if fdownstream_idx == len(tss['fwd_genes']):
                   fdownstream_idx = fdownstream_idx -1
                   rup_ofs = 0
               
               #get gene indexes in the chromosomal gene dicts
               #for upstream and downstream
               fdown_gene = tss['fwd_genes'][fdownstream_idx]
               fup_gene = tss['fwd_genes'][fdownstream_idx +fup_ofs]\
                   if fdownstream_idx > 0 \
                   else fdown_gene
               
               rdown_gene = tss['rev_genes'][rdownstream_idx]
               rup_gene = tss['rev_genes'][rdownstream_idx +rup_ofs]\
                   if rdownstream_idx < len(tss['rev_genes']) - 1\
                   else rdown_gene
             
               bind_infos[tfname][ekey].append({
                       'rup_gene':rup_gene,
                       'rdown_gene':rdown_gene,
                       'fup_gene':fup_gene,
                       'fdown_gene':fdown_gene,
                       'chr':crkey,
                       'start':start,
                       'end':end,
                       'mean':(start+end)/2,
                       'score':e['score']                  
       
                })                 
       
       return bind_infos
    return mem.getOrSet(set_tf_chips, **mem.rc(kwargs))

def get_assay_info(**kwargs):
    kwargs['atype'] = kwargs.get('atype', default_atype)
    def set_assay_info(**kwargs):
       #peaks and keys 
       atype = kwargs['atype']
       if atype == 'tfchip':
           assay_peaks = tf_chip_peaks(**mem.sr(kwargs))
       elif atype=='wormtile':
           assay_peaks = tiling_peaks(**mem.sr(kwargs))
       else: raise Exception()
       pkeys = assay_peaks.keys()
       
       keyvalfuns = {'tfchip':lambda x: x[:x.index(':')],
                     'wormtile':lambda x:re.compile('Tissue=([^\(]*)')\
                         .search(x).group(1)}
       #index and uniquify tfs
       pkeyvals = dict([(k, keyvalfuns[atype](k)) for k in pkeys])
       tfnames = set(pkeyvals.values())
       tfkeys= dict([(name, [k for k,v in pkeyvals.iteritems()
                             if v == name])
                     for name in tfnames])       

       chrmap = dict(zip(['I','II', 'III','IV', 'V', 'X'],
                         chromosome_names()))
       
       all_tss = get_tss()
       all_genes = parse_genes()
       
       bind_infos = {}
       #loop through tfs, assays, and finally peaks
       for tfname in tfnames:
           print tfname
           bind_infos[tfname] = {}
           for ekey in tfkeys[tfname]:
             print 'n_exps = {0}'.format(np.sum([len(v) 
                                              for v in bind_infos.values()]))
             exp = assay_peaks[ekey]
             bind_infos[tfname][ekey] = []
             for e in exp:
               start = e['start']
               end = e['end']
               c = e['chr']

               #if we have non chromosomal DNA, skip it!
               if c in chrmap.keys():
                   crkey = chrmap[c]
               else: continue
               genes = all_genes[crkey]
               tss = all_tss[crkey]
               
               #grab indexes in the list of upstream and ds genes
               fdownstream_idx =searchsorted(tss['fwd_tss'],(start+end) /2)
               rdownstream_idx =searchsorted(tss['rev_tss'],(start+end) /2, 'right')-1
               
               #handle the case where the factor hits past the 
               #last gene
               fup_ofs = -1
               if fdownstream_idx == len(tss['fwd_genes']):
                   fdownstream_idx = fdownstream_idx -1
                   fup_ofs = 0
               
               rup_ofs = 1
               if fdownstream_idx == len(tss['fwd_genes']):
                   fdownstream_idx = fdownstream_idx -1
                   rup_ofs = 0
               
               #get gene indexes in the chromosomal gene dicts
               #for upstream and downstream
               fdown_gene = tss['fwd_genes'][fdownstream_idx]
               fup_gene = tss['fwd_genes'][fdownstream_idx +fup_ofs]\
                   if fdownstream_idx > 0 \
                   else fdown_gene
               
               rdown_gene = tss['rev_genes'][rdownstream_idx]
               rup_gene = tss['rev_genes'][rdownstream_idx +rup_ofs]\
                   if rdownstream_idx < len(tss['rev_genes']) - 1\
                   else rdown_gene
             
               bind_infos[tfname][ekey].append({
                       'rup_gene':rup_gene,
                       'rdown_gene':rdown_gene,
                       'fup_gene':fup_gene,
                       'fdown_gene':fdown_gene,
                       'chr':crkey,
                       'start':start,
                       'end':end,
                       'mean':(start+end)/2,
                       'score':e['score']                  
       
                })                 
       
       return bind_infos
    return mem.getOrSet(set_assay_info, 
                        **mem.rc(kwargs ,
                                 name = kwargs['atype']))


def get_assay_gprops(**kwargs):
    kwargs['atype'] = kwargs.get('atype', default_atype)
    def set_assay_gprops(**kwargs):
       chips = get_assay_info(**mem.sr(kwargs))
       genes = parse_genes()
       tf_stats = {}
       for k, v in chips.iteritems():
           tf_stats[k] = {}
           for k2,v2 in v.iteritems():
               print 'n_exps = {0}'.\
                   format(np.sum([len(v) 
                                  for v in tf_stats.values()]))

               tf_stats[k][k2] = {}
               cs = [e['chr'] for e in v2]
               f_gups=[genes[cs[i]][e['fup_gene']] for i,e in enumerate(v2)]
               f_gdowns=[genes[cs[i]][e['fdown_gene']] for i,e in enumerate(v2)]
               
               fup_deltas = [e['mean'] - f_gups[i].location.start.position
                             for i,e in enumerate(v2)]
               fdown_deltas=[e['mean'] - f_gdowns[i].location.start.position
                             for i,e in enumerate(v2)]
               
               
               r_gups=[genes[cs[i]][e['rup_gene']] for i,e in enumerate(v2)]
               r_gdowns=[genes[cs[i]][e['rdown_gene']] for i,e in enumerate(v2)]
               
               rup_deltas = [e['mean'] - r_gups[i].location.end.position
                             for i,e in enumerate(v2)]
               rdown_deltas=[e['mean'] - r_gdowns[i].location.end.position
                             for i,e in enumerate(v2)]
               
               deltas = array([fdown_deltas,
                               fup_deltas,
                               rdown_deltas,
                               rup_deltas]).T
               
               closest = argmin(np.abs(deltas),1)
               csrt = argsort(np.abs(deltas),1)
               
               primaries = []
               secondaries = []
               for i,c in enumerate(csrt):
                   for j,e in enumerate(c[:2]):
                       arr = primaries if j == 0 else secondaries
                       d = {}
                       if e == 0 : d['gene'] = f_gdowns[i]
                       elif e==1 : d['gene'] = f_gups[i]
                       elif e==2 : d['gene'] = r_gdowns[i]
                       elif e==3 : d['gene'] = r_gups[i]
                       
                       d['dist'] = deltas[i,e] * (-1 if e < 2 else 1)
                       arr.append(d)
       
               tf_stats[k][k2]['primaries'] = primaries
               tf_stats[k][k2]['secondaries'] = secondaries
       return tf_stats
    return mem.getOrSet(set_assay_gprops, **mem.rc(kwargs,
                                                    name = kwargs['atype'],
                                                    hardcopy = True))

def get_simple_description(**kwargs):
    kwargs['atype'] = kwargs.get('atype', default_atype)
    def set_simple_description(**kwargs):
        props = get_assay_gprops(**mem.sr(kwargs))
        chips = get_assay_info(**mem.sr(kwargs))
        simple = {}
        for tf, assays in props.iteritems():
            simple[tf] = {}
            assay_keys = assays.keys()

            idfun = lambda g: g.qualifiers['db_xref'][1][9:] 

            simple[tf]['gnames'] = list(it.chain(*[ 
                        [idfun(e['gene']) for e in assays[k]['primaries']]
                        for k in assay_keys
                        ]))
            simple[tf]['genes'] = list(it.chain(*[ 
                        [e['gene'] for e in assays[k]['primaries']]
                        for k in assay_keys
                        ]))
            simple[tf]['dists'] = array(list(it.chain(*[
                        [e['dist'] for e in assays[k]['primaries']]
                        for k in assay_keys
                        ])))
            simple[tf]['scores'] = array(list(it.chain(*[
                        [e['score'] for e in chips[tf][k]]
                        for k in assay_keys
                        ])))
        return simple
    return mem.getOrSet(set_simple_description,
                        **mem.rc(kwargs,
                                 name = kwargs['atype']))
                                                 
    
def get_simple_thr(**kwargs):
    kwargs['dthr'] = kwargs.get('dthr',None)
    kwargs['sthr'] = kwargs.get('sthr',None)
    kwargs['dsign'] = kwargs.get('dsign',None)
    kwargs['atype'] = kwargs.get('atype', default_atype)

    def set_simple_thr(**kwargs):
        dthr = kwargs['dthr']
        dsign = kwargs['dsign']
        sthr = kwargs['sthr']
        simple = get_simple_description(**mem.sr(kwargs))
        out = {}
        for k,v in simple.iteritems():
            criteria = ones(len(v['scores']))
            if sthr != None:
                criteria *= less(v['scores'],sthr)
            if dsign != None:
                criteria *= greater(v['dists']*dsign, 0)
            if dthr != None:
                criteria *= less(abs(v['dists']),dthr)
                            
            allowed = nonzero(criteria)[0]
            out[k] = {'genes':[v['genes'][i] for i in allowed],
                      'gnames':[v['gnames'][i] for i in allowed],
                      'dists':v['dists'][allowed],
                      'scores':v['scores'][allowed]}
        return out

    #names = {'wormtile':'sthr_{0}'.format(kwargs['sthr']),
    #         'tfchip':'dthr_{0}_sthr_{1}'.\
    #             format(kwargs['dthr'],
    #                    kwargs['sthr'])
    #         }
    
    tkwargs = dict([(k,kwargs[k]) for k in ['dthr','sthr', 'dsign']])
    name = '{0}:'.format(kwargs['atype']) + \
        '_'.join(it.chain(*sorted([(str(k),str(v)) for k,v in tkwargs.iteritems()],
                                  key = lambda x: x[0])))
    
                                  
    return mem.getOrSet(set_simple_thr,  
                        **mem.rc(kwargs,
                                 name =name))
                                                 
    




def tiling_peaks(**kwargs):
    def set_tiling_peaks(**kwargs):
        root = cfg.dataPath('modencode/wormtile/computed-peaks_gff3')
        files = [os.path.join(root, f) for f in os.listdir(root)]
        out = {}
        for f in files:
            if f[-2:] != 'gz': continue
            fopen= gzip.open(f)
            data = [l for l in fopen.readlines() if not l[0] == '#']
            
            out[os.path.basename(f)] = \
                [dict(zip(['chr', 'meth', 'type',
                  'start','end','score',
                  'blank','blank2','annotations' ], 
                         l.strip().split('\t')))
                 for l in data]
        
        for k,v in out.iteritems():
            for d in v:
                d['start'] = int(d['start'])
                d['end'] = int(d['end'])
                d['score'] = float(d['score'])
        return out
    return mem.getOrSet(set_tiling_peaks, 
                        **mem.rc(kwargs,
                                 hardcopy = True,
                                 name = 'default'))


def chromosome_names():
    return  ['CHR_I', 'CHR_II', 'CHR_III'
             ,'CHR_IV','CHR_V','CHR_X']
def chromosome_offsets(**kwargs):
    def set_chromosome_offsets(**kwargs):

      lens =[]
      names = chromosome_names()
      for name in names:
         root = cfg.dataPath('/data/genomes/Caenorhabditis_elegans')
         fdir = os.path.join(root,name)
         for r, d, files in os.walk(fdir):
             for f in files:
                 if '.gb' in f:
                     fopen = open(os.path.join(r,f))
                     break
      
         gb = list(sio.parse(fopen, 'genbank'))[0]
         fopen.close()
         lens.append( gb.features[0].location.end.position)

      offsets = {}
      cur_ofs = 0
      for i, l in enumerate(lens):
        offsets[names[i]] = cur_ofs
        cur_ofs += l 
      return offsets
    return mem.getOrSet(set_chromosome_offsets, 
                        **mem.rc(kwargs, hardcopy = True))
    
def parse_genes(**kwargs):

    def set_genes(**kwargs):

      lens =[]
      all_genes = {}
      names = chromosome_names()
      for name in names:
         root = cfg.dataPath('/data/genomes/Caenorhabditis_elegans')
         fdir = os.path.join(root,name)
         for r, d, files in os.walk(fdir):
             for f in files:
                 if '.gb' in f:
                     fopen = open(os.path.join(r,f))
                     break
      

         gb = list(sio.parse(fopen, 'genbank'))[0]
         genes = [f for f in gb.features if f.type== 'gene']
         all_genes[name] = genes
         fopen.close()
      return all_genes
    return mem.getOrSet(set_genes, 
                        **mem.rc(kwargs,hardcopy = True))

def cr_locii(chrname = 'CHR_I', **kwargs):
    ng =20;
    
    root = cfg.dataPath('/data/genomes/Caenorhabditis_elegans')
    fdir = os.path.join(root,chrname)
    for r, d, files in os.walk(fdir):
        for f in files:
            if '.gb' in f:
                fopen = open(os.path.join(r,f))
                break
    gb = list(sio.parse(fopen, 'genbank'))[0]

    chrinfos= parse_genes();
    for i in range(ng):
        g = chrinfos[chrname][i]
        gid = g.qualifiers['db_xref'][1][9:] 
        print gid
        neighbor_range =[max(0, g.location.start.position-5000),
                         min(len(gb.seq)-1, g.location.end.position+5000)]
        feats_matching =[f for f in gb.features
                         if f.location.start.position > neighbor_range[0]
                         and f.location.end.position < neighbor_range[1]]
        
        seq = gb.seq[neighbor_range[0]:neighbor_range[1]]
    
        import Bio.SeqRecord as sr
        import Bio.Seq as sq
        rec = sr.SeqRecord(sq.Seq(seq.tostring(),seq.alphabet),
                           name = '{0}_g{1}'.format(chrname,i),
                           description='genomic neighborhood of {0}'\
                               .format(gid),
                           annotations = {'organism':'worm'})

        import Bio.SeqFeature as sf
        import copy
        for f in feats_matching:
            
            fnew = copy.deepcopy(f)
            fnew.location.start.position -= neighbor_range[0]
            fnew.location.end.position -= neighbor_range[0]
            rec.features.append(fnew)
        root = cfg.dataPath('modencode/worm_gene_extracts/{0}'.format(chrname))
        if not os.path.isdir(root): os.mkdir(root)
        fname = os.path.join(root,'{0}.gbk'.format(gid))
        fopen = open(fname, 'w')
        sio.write(rec, fopen, 'genbank')
        fopen.close()
        
        


def gene_info(**kwargs):
    def set_gene_info(**kwargs):
        crofs = chromosome_offsets()
        gp = parse_genes()
        genes = {}
        idfun = lambda g: g.qualifiers['db_xref'][1][9:] 
            
        for cname, elts in gp.iteritems():
            crof = crofs[cname]
            genes.update(dict([(idfun(e), {'start':e.location.start.position,
                                           'end':e.location.end.position,
                                           'genomestart':e.location.start.position+crof,
                                           'genomeend':e.location.end.position+crof,
                                           'tss':e.location.start.position  \
                                               if e.strand == 1 \
                                               else e.location.end.position})
                               for e in elts]))
        return genes
    return mem.getOrSet(set_gene_info,**mem.rc(kwargs,
                                               name = 'default'))                                         
                                         
                                           

def get_tss(**kwargs):
    def load_tss(**kwargs):
       cnames = chromosome_names()
       genes = parse_genes()
       
       out = {}
       for name in cnames:
           crgenes= genes[name]
           gstrands =array([g.strand for g in crgenes])
           fwd = nonzero(greater(gstrands, 0))[0]
           rev = nonzero(less(gstrands,0))[0] 
           
           gstarts = array([g.location.start.position for g in crgenes])
           gends =array( [g.location.end.position for g in crgenes])
       
           fstarts = gstarts[fwd]
           rends = gends[rev]
           
           fstart_sorted = sorted([(fwd[i], s)
                                   for i, s in enumerate(fstarts)],
                                  key = lambda x: x[1])
           fend_sorted = sorted([(rev[i], r)
                                 for i, r in enumerate(rends)],
                                key = lambda x: x[1])
           
           out[name] = {'fwd_genes': [e[0] for e in fstart_sorted],
                        'fwd_tss':[e[1] for e in fstart_sorted],
                        'rev_genes': [e[0] for e in fend_sorted],
                        'rev_tss':[e[1] for e in fend_sorted]}
       
           #note that gstarts begin sorted
           #gends on the other hand... do not.
       return out
    return mem.getOrSet(load_tss, **mem.rc(kwargs,hardcopy =True))
