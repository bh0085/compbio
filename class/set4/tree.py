import matplotlib.pyplot as plt
from numpy import *
import matplotlib.text
import numpy as np
import upmg

v = True

seed_seqs = ['AACCGG',
             'ACTCAG',
             'GTCCTT',
             'GGTTCG']
seq_translate = {'A':0,'G':1,'T':2,'C':3}
costs =array( [[0,1,2,2],
        [1,0,2,2],
        [2,2,0,1],
        [2,2,1,0]
        ])

class Tree2():
    def __init__(self, parents, children, labels = None, weights = None):
        self.parents = parents
        self.children = children
        self.labels = labels
        self.weights = weights
        print weights
        self.nodexy = None
        self.node_xy = [[] for i in range(len(parents))]
        if v: print 'canonicalizing tree on init'
        self.canonicalize()
        self.make_xys()
        self.cxnlabels = None
        self.biglabel = ''


    def canonicalize(self):
        p2 = list(self.parents)
        c2 = list(self.children)
        #canonicalizes the list of node children.
        print p2
        #root the tree at the last node in the list:
        c0 = c2[-1]
        
        if v: print "canonicalizing under the assumption that nodes are"
        if v: print "indexed with children appearing before parents"

        n = len(p2)
        min_child = arange(n)
        for i in range(n):
            min_child[p2[i]] = min([min_child[p2[i]],min_child[i]])
            
    
        print c2
        for i in range(n):
            if not c2[i]: continue
            if min_child[c2[i][0]] > min_child[c2[i][1]]:
                c2[i].reverse()
        self.children = c2
        if v: print "done canonicalizing"
        

    def make_xys(self):
        leaves = nonzero( map(lambda x: x == [], self.children))[0]
        n = len(self.children)
        sorted_leaves = []
        self.subtree_leaves(n-1,sorted_leaves)
        nl = len(leaves)
        
        #generate a list of locations and spool them into
        #leaf coordinates ordered in the canonical fashion
        xs = 10*cos((arange(float(nl))-0)/nl * pi * 2)
        ys = 10*sin((arange(float(nl))-0)/nl * pi * 2)
        
        self.lx, self.ly = [], []
        e_idx = -1
        for i in range(nl):
            try: 
                this_i = sorted_leaves.index(i)
            except Exception, e:
                this_i = e_idx
                e_idx -= 1


            self.lx.append(  xs[this_i] )
            self.ly.append(  ys[this_i] )

        roots = nonzero(equal(self.parents,-1))[0]
        for r in roots:
            xy_root = self.get_xy(r)

        print self.node_xy
        
    def pull_xys(self, pull = .4, leaves = True):
        n = len(self.children)
        for i in range(n)[::-1]:
            xyp = self.node_xy[i]
            for c in self.children[i]:
                if not leaves and self.children[c] == []: continue
                self.node_xy[c] += (xyp - self.node_xy[c])*pull
                

    
    def scale_branches(self):
        print 'scaling branches with a pad of 1'
        pad = 1
        roots = nonzero(equal(self.parents,-1))[0]
        for r in roots:
            self._scale_branches(r,[0.0,0.0], pad = 1)

    def _scale_branches(self,idx, shiftxy,pad = 0):
        """shiftxy is the amount that the parent of this node has already been shifted"""

        childs = self.children[idx]
        weights = self.weights[idx]

        newxy =  self.node_xy[idx] + shiftxy
        for i in range(len(childs)):

            
            c = childs[i]
            w = weights[i]
            cxy = self.node_xy[c]
            delta = cxy - newxy
            mag = np.sqrt(np.sum(np.power(delta,2)))
            cshift = (newxy + (delta/mag *(w+ pad))) - cxy 
            self._scale_branches(c, cshift,pad = pad)

        self.node_xy[idx] = newxy
        

        # + delta / mag * w 

    def subtree_leaves(self,idx,sorted_leaves):
        """return a list of subtree leaves in the order of their occurence"""
        #... using a recursive alg...
        if not self.children[idx]:
            sorted_leaves.append(idx)
            return
        else:
            self.subtree_leaves(self.children[idx][0],sorted_leaves)
            self.subtree_leaves(self.children[idx][1],sorted_leaves)
            return

    def get_xy(self,idx):
        if not self.node_xy[idx]:
            self.compute_xy(idx)
        return self.node_xy[idx]
        
    def compute_xy(self, idx):
        if self.children[idx] == []:
            self.node_xy[idx] =array( [self.lx[idx],self.ly[idx]])
        else:
            self.node_xy[idx] = np.mean([self.get_xy(self.children[idx][0]),
                                    self.get_xy(self.children[idx][1])], 0)
        

def run_upgma(max_itr = -1):
    max_n = 4
    f = plt.figure(0)
    f.clear()
    f.set_facecolor('.8')

    for i in range(max_n):
        ax = f.add_subplot('22'+str(i+1))
        u,dists  = upmg.run_test_upgma(max_itr = i+2)
        t = u.tree
        t.biglabel = str(dists)
        t.pull_xys(.25, leaves = False)
        draw(t, plot_weights = False, ax = ax)
        
    ax.annotate('(b) UPGMA with distance matrices',[0,0],xytext = [.05,.9],family = 'serif',textcoords = 'figure fraction',size = 'xx-large')
        
        

def run_nj(max_itr = -1):

    max_n = 4
    f = plt.figure(0)
    f.clear()
    f.set_facecolor('.8')

    
    ax = f.add_subplot('11'+str(1))
    u,dists = upmg.run_test_nj()
    t = u.tree
    t.scale_branches
    t.pull_xys(.4, leaves = False)

    draw(t, plot_weights = True, ax = ax)
    ax.annotate('(c) The NJ Algorithm Get Lengths Right',[0,0],xytext = [.05,.9],family = 'serif',textcoords = 'figure fraction',size = 'xx-large')


def run():

    coals = [[[0,1],[0,1],[0,1]],
             [[0,2],[0,1],[0,1]],
             [[0,3],[0,1],[0,1]]]


    f = plt.figure(0)
    f.clear()
    f.set_facecolor('.8')




    for i in range(len(coals)):
        c = coals[i]
        ax = f.add_subplot('13'+ str(i))
        t = build(nodes = range(4),coals = c)
        t.pull_xys(leaves = False)
        s = get_costs(t)
        draw(t,plot_weights = False,ax = ax)
    ax.annotate('(a) The best pairing is the third one',[0,0],xytext = [.05,.9],family = 'serif',textcoords = 'figure fraction',size = 'xx-large')

def get_coals():
    return [[0,1],[0,1],[0,1]]
def get_costs(tree):
    c = tree.children
    n = len(c)
    leaves = nonzero(map(lambda x : x == [] and True or False , c))[0]
    int_nodes = nonzero(map(lambda x: x!= [] and True or False,c))[0]
    
    seen = zeros(n,int)
    seen[leaves] = 3
    
    compute = []
    print 'assuming fixed seq_length'
    seq_len = len(tree.labels[0])
    print seq_len
    seq_costs = zeros((n,seq_len,4),float)

    
    for l in leaves:
        seen[tree.parents[l]]+=1
        s = map(lambda x: seq_translate[x],tree.labels[l])
        for i in range(seq_len):
            seq_costs[l][i][:] = -1
            seq_costs[l][i][s[i]] = 0
    for w in nonzero( equal(seen,2))[0]:
        seen[w] = 3
        compute.append(w)

    #print seq_costs
    while len(compute):

        for c in compute:
            ch = tree.children[c]
            for i in range(seq_len):
                for j in range(4):
                    seq_costs[ch[0],:]
                    idx0 = nonzero(greater(seq_costs[ch[0],i],-1))[0]
                    idx1 = nonzero(greater(seq_costs[ch[1],i],-1))[0]
                    c0_costs = np.min(seq_costs[ch[0],i,idx0] + costs[j,idx0])
                    c1_costs = np.min(seq_costs[ch[1],i,idx1] + costs[j,idx1])
                    seq_costs[c,i,j]=c0_costs + c1_costs
            compute.remove(c)
            seen[tree.parents[c]] +=1
        for w in nonzero( equal(seen,2))[0]:
            seen[w] = 3
            compute.append(w)
    
    labels = []
    max_c = -1
    for i in range(n):
        label = []
        for j in range(seq_len):
            p_costs = np.min(seq_costs[i,j,:])
            label.append(p_costs)
        if np.sum(label) > max_c: max_c = np.sum(label)
        labels.append(label)
    for i in range(len(tree.labels)):
        if tree.labels[i] == '':
            tree.labels[i] = labels[i]
    tree.biglabel = 'Cost: ' +str(max_c)
def build(nodes = range(4), coals = [[0,1],[0,1],[0,1]]):
    #use a fairly silly notation for coalescence
    #store the current cluster of each leaf so that
    #I can merge leaves...
    currents = range(4)
    nlinks = []
    seqs = []
    for x in nodes:
        nlinks.append([])
        seqs.append(seed_seqs[x])

    T=nodes
    L=list(T)
    parents = [-1] * len(T)
    children = [[]] * len(T)
    weights = list(children)
    for c in coals:
        parents.append(-1)
        children.append([])
        weights.append([1,1])
        seqs.append('')
        nidx = len(T)
        T.append(nidx)

        chosen= [L[c[0]],L[c[1]]]
        for n in chosen:
            L.remove(n)
            parents[n] = nidx
            children[nidx].append(n)
        print nidx
        print T
        L.append(nidx)
        

    t = Tree2(parents, children, labels = seqs,weights = weights)

    return t



def merge_seqs(s1,s2):
    return s1

def node_radius(tree,idx):
    if tree.children[idx] == []:
        return 20
    else:
        return 5

def draw(tree,plot_weights = True, ax = None):

    xs = []
    ys = []
    rs = []
    ls = []
    l_ofs = []
    import matplotlib.patches as patches
    cxns = []
    n = len(tree.parents)
    cxweights = []

    for i in range(n)[::-1]:
        children = tree.children[i]
        

        if tree.cxnlabels:
            labels = tree.cxnlabels[i]
        if tree.weights:
            labels = map(lambda x: str(around(x,3)),tree.weights[i])
        else:
            labels = ['' for i in range(len(children))]
            

        x0, y0 = tree.node_xy[i]
        xs.append(x0)
        ys.append(y0)
        
        ls.append(tree.labels[i])

        r = node_radius(tree,i)
        rs.append(pow(r,2)*pi)
        l_ofs.append(r+3)         
        for i in range(len(children)):
            w = labels[i]
            c = children[i]
            x1,y1 = tree.node_xy[c]
            delta = array([x0,y0]) - array([x1,y1])
            cxweights.append({'xy':np.mean(array([[x0,y0],[x1,y1]]),0),
                              'wt':w,
                              'unit':delta / sqrt(np.sum(power(delta,2)))})

            p = patches.ConnectionPatch([x0,y0],[x1,y1],'data','data',
                                        arrowstyle = 'wedge,tail_width=1.5',
                                        shrinkA = r + 2,
                                        shrinkB = node_radius(tree,c) + 3,
                                        edgecolor = 'black',
                                        facecolor = 'white',
                                        alpha = .5)
            cxns.append(p) 
        
    if not ax:
        f = plt.figure(0)
        f.clear()
        ax = f.add_subplot('111', frameon = False)
        
    plt.axis('off')


    for x in cxns:
        ax.add_patch(x)


    ax.scatter(xs,ys,rs,'white')
    for i in range(len(xs)):
        ax.annotate(ls[i],[xs[i],ys[i]],
                    xytext = [l_ofs[i],l_ofs[i]],
                    textcoords = 'offset points',
                    bbox=dict(boxstyle="round4,pad=.5", fc="0.99",alpha = .5),
                    color = (lambda x: x and 'black' or 'black')(plot_weights),
                    size = 'small')

    if plot_weights:
        for w in cxweights:
            print sqrt(sum(power(w['unit'],2)))
            xytext = 50*array([w['unit'][1],-1*w['unit'][0]])
            print  sqrt(sum(power(xytext,2)))
            ax.annotate(w['wt'],w['xy'], 
                        xytext = xytext,
                        textcoords = 'offset points',
                        arrowprops = {'width':10,'alpha':1,'shrink':.1},
                        color = 'black',
                        size = 'small')
        
    if tree.biglabel:
        ax.annotate(tree.biglabel, [0,0],
                    family = 'serif',size = 'x-large',
                    xycoords='axes points')
