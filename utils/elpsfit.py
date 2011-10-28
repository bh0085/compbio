import numpy as np
from numpy import *
import itertools as it

def cluster_ellipses(coords, tis):
        
    tis_means = array([np.mean([coords[:,idx[0]] for idx in g],0)  
                 for k, g in it.groupby(enumerate(tis),
                                        key = lambda x: x[1])])
    tis_vars = [ np.var([coords[:,idx[0]] for idx in g],0)
                 for k, g in it.groupby(enumerate(tis),
                                        key = lambda x: x[1])]
    centered =array( [coords[:,i] - tis_means[tis[i],:] for i in range(len(tis))])

    tis_cov = [sum([ dot(centered[idx[0],:][:,newaxis],
                         centered[idx[0],:][newaxis,:]) for idx in g],0) 
               for k, g in it.groupby(enumerate(tis),
                                      key = lambda x: x[1])]   
    tis_eigs = [linalg.eig(x) for x in tis_cov]

    import matplotlib.patches as patch

    elpses = []
    chars = []
    for i, elt  in enumerate(zip(tis_means, tis_vars, tis_eigs)):
        m, v ,e= elt
        vals, vecs = e
        
        mean = m[:2] 
        vars = array(vals[:2])
        vecs  =  array(vecs[:2])
        
        
        vars = sqrt(vars)
        chars.append({'mean':mean,
                      'axis_lengths':vars,
                      'axis_vectors':vecs})
        el = patch.Ellipse(mean, *(vars/10),\
                               angle = 1.*180  / pi *\
                               arctan(vecs[0][1]/vecs[0][0]) )
        elpses.append(el)

    return elpses, chars
