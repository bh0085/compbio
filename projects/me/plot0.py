import cb.config as cfg
import numpy as np
from numpy import *
import os
import cb.utils.memo as mem
import cb.utils.graphs.draw as gd

import networkx as nx

import utils as wu
import io as io
import parse as wp
import inference as wi
import cb.utils.plots as myplots

def peak_distance_histogram(**kwargs):    

    atype = kwargs.get('atype', wp. default_atype)
    chips = wp.get_assay_gprops(**mem.rc(kwargs))
    chiplist = chips.values()
    chipkeys = chips.keys()
    xs = []
    ys = []
    

    sec_spread = np.max([
            np.max([ 
                    np.max(np.abs([e['dist'] for e in v2['secondaries']]))
                    for v2 in v.values()])
            for v in chips.values()
            ])

    hist_spread = 10000
    bin_wid = 200
    bin_mids = arange(-1* hist_spread, 1*hist_spread,bin_wid)
    bin_starts = bin_mids - bin_wid/2
    nb = len(bin_starts)

    prim_hists = zeros((len(chips), len(bin_starts)))
    sec_hists = zeros((len(chips), len(bin_starts)))
    for i,e in enumerate(chiplist):
        for k,v in e.iteritems():
            pbins = array([e2['dist']/bin_wid for e2 in v['primaries']],int)
            sbins = array([e2['dist']/bin_wid for e2 in v['secondaries']],int)
            pbins += nb /2
            sbins += nb /2

            sbins[less(sbins,0)] = 0
            sbins[greater(sbins,nb-1)] = nb-1
            
            pbins[less(pbins,0)] = 0
            pbins[greater(pbins,nb-1)] = nb-1
            for b in pbins: prim_hists[i][b]+=1
            for b in sbins: sec_hists[i][b]+=1
            
    f= myplots.fignum(1,(8,6))
    ax = f.add_subplot(111)
    ax.set_title('chip peak distances to primary/sec tss for {0}'.format(atype))
    for p in prim_hists:
        ax.plot(bin_mids,p, color = 'green')
    for s in sec_hists:
        ax.plot(bin_mids,s, color = 'red')
    f.savefig(myplots.figpath('chip_distance_hists_for{0}.pdf'.format(atype)))

    

def peak_thr_histograms(**kwargs):
    '''histograms of score and distance'''
    dthr = 1500
    sthr = 1e-2
    dsign = -1
    simple = wp.get_simple_thr(**mem.rc(kwargs,
                                        dthr = dthr,
                                        sthr = sthr,
                                        dsign = dsign
                                        )
                               )

                                 
                                 
    min_score = -1
    max_score = -1
    for k,v in simple.iteritems():
        smax = np.max(v['scores'])
        if max_score == -1 or smax > max_score:
            max_score = smax
        smin = np.min(v['scores'])
        if min_score == -1 or smin < min_score:
            min_score = smin

    lrange = [int(floor(log10(min_score))), int(ceil(log10(max_score)))]
    sbin_mids = range(lrange[0],lrange[1]+1)
    nsb = len(sbin_mids)
    sbins = zeros((nsb))

    dbin_size = 50
    dbin_mids = range(-dthr, dthr, dbin_size)
    ndb = len(dbin_mids)
    dbins = zeros(( ndb))

    for k,v in simple.iteritems():
        for d in v['dists']:
            dbins[int(d + dthr)/dbin_size] += 1
        for s in v['scores']:
            sbins[int(log10(s) - lrange[0])] += 1
    
    f= myplots.fignum(1,(10,6))
    ax = f.add_subplot(121)
    ax.set_ylabel('log 10 counts')
    ax.set_xlabel('distance')
    ax.set_title('simplified tss distances (d<{0})'.format(dthr))
    ax.plot(dbin_mids,log10(dbins), color = 'black')
    f.savefig(myplots.figpath('chip_simple_distance.pdf'))

    ax = f.add_subplot(122)
    ax.set_title('simplified tss scores (s<{0:2.2})'.format(sthr))
    ax.set_ylabel('log10 counts')
    ax.set_xlabel('log10 peak score')
    ax.plot(sbin_mids,log10(sbins), color = 'black')
    f.savefig(myplots.figpath('chip_simple_scores.pdf'))


def plot_easy_inference():
    dg = io.getGraph()
    pos = gd.getpos(dg)
    
    f = myplots.fignum(4, (8,8))
    ax = f.add_subplot(111)
    ax.set_title('putative worm chip network')
    gd.easy_draw(dg, pos)

    f.savefig(myplots.figpath('worm_chip_graph.pdf'))


def plot_coexpression():
    genes, gene_info = wi.get_coexpression()
    
