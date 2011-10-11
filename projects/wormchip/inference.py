import parse as wp
from numpy import *
import numpy as np
import itertools as it

import utils as wu
import cb.utils.memo as mem

def get_easy0(**kwargs):
    '''an easy inference using the arbitrary distance cutoff of 3000 bases
    and grabbing the n highest scoring edges globally'''
    def set_easy0(**kwargs):
       simple = wp.tf_chip_simple_thr()
       score_soft_cut = -136
       score_hard_cut = -90
       sids = wu.symbol_ids()
       prop_tuples = []
       for tf ,props in simple.iteritems():
           #for now, remove tfs that are not mappable
           if not tf in sids.keys(): continue
           ssrt = argsort(props['scores'])
           lscores = log10(props['scores'][ssrt])
       
           easy = nonzero(less(lscores , score_soft_cut))[0]
           medium = nonzero(greater(lscores, score_soft_cut)*\
                                less(lscores, score_hard_cut))[0]
           #generous_edges = concatenate([easy,medium])
           prop_tuples.append([(tf, props['genes'][g] ,
                                -(score_hard_cut - lscores[g]) \
                                    /(score_soft_cut - score_hard_cut)) 
                               for g in medium])
           prop_tuples.append([(tf, props['genes'][g],1) for g in easy])
           
       edgelist = array(list(it.chain(*prop_tuples)))
       edges = [( sids[e[0]] ,e[1].qualifiers['db_xref'][1][9:],e[2] )
              for e in edgelist]
       return edges

    return mem.getOrSet(set_easy0, **mem.rc(kwargs))
