import matplotlib.pyplot as plt
from numpy import *
import matplotlib.text
import numpy as np

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

class MyTree():
    def __init__(self,nx, ny, nlinks,parents, labels = []):
        self.xs = nx
        self.ys = ny
        self.children = nlinks
        self.parents = parents
        self.labels = labels

def run():
    t = build()
    s = get_costs(t)
    print t.labels
    draw(t)

def get_costs(tree):
    c = tree.children
    n = len(c)
    leaves = nonzero(map(lambda x : x == [] and True or False , c))[0]
    int_nodes = nonzero(map(lambda x: x!= [] and True or False,c))[0]
    
    seen = zeros(n,int)
    seen[leaves] = 3
    
    compute = []
    print 'fixing seq length = 6'
    seq_len = 6
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
    for i in range(n):
        label = []
        for j in range(seq_len):
            label.append(np.min(seq_costs[i,j,:]))
        labels.append(label)
    for i in range(len(tree.labels)):
                                   
        if tree.labels[i] == '':
            tree.labels[i] = labels[i]

def get_coals(x,y):
    return [[3,1],[2,0],[0,3]]

def build(nodes = range(4)):


    xleaves = array([0,0,1,1],float)
    yleaves = array([0,1,0,1],float)
    
    #use a fairly silly notation for coalescence
    #store the current cluster of each leaf so that
    #I can merge leaves...

    coals = get_coals(xleaves, yleaves)
    currents = range(4)
    
    nlinks = []
    nx = []
    ny = []
    seqs = []
    for x in nodes:
        nlinks.append([])
        nx.append(xleaves[x])
        ny.append(yleaves[x])
        seqs.append(seed_seqs[x])

    for i in coals:
        i_last = i
        #search forward for node current positions in the list
        while 1:
            i_new= currents[i_last[0]],currents[i_last[1]]
            if i_new == i_last: break
            else: i_last = i_new
        idxs = i_new

        x = mean(array([nx[idxs[0]],nx[idxs[1]]],float))
        y = mean(array([ny[idxs[0]],ny[idxs[1]]],float))
        
        seqs.append('')

        cur_node = len(nlinks)

        #link the list backward
        nlinks.append(idxs)
        nx.append(x)
        ny.append(y)


        #link the list forward
        currents.append(cur_node)

        currents[idxs[0]] = cur_node
        currents[idxs[1]] = cur_node

    t = MyTree(nx,ny,nlinks,currents, labels = seqs)
    return t

def merge_seqs(s1,s2):
    return s1

def draw(tree):
    #PARAMS
    nudge = .25

    xs = []
    ys = []
    rs = []
    ls = []
    l_ofs = []
    
    import matplotlib.patches as patches
    cxns = []

    node_xs = array(tree.xs)
    node_ys = array(tree.ys)

    #Pull Children towards parents
    n = len(node_xs)
    for i in range(n)[::-1]:
        p = i
        children = tree.children[i]
        for c in children: 
            node_xs[c] +=  nudge *( node_xs[p] - node_xs[c])
            node_ys[c] +=  nudge *( node_ys[p] - node_ys[c])


    for i in range(n)[::-1]:
        children = tree.children[i]
        p = tree.parents[i]

        x0 = node_xs[i]
        y0 = node_ys[i]
        xs.append(x0)
        ys.append(y0)
        
        ls.append(tree.labels[i])


        if children:
            r = 5
            rs.append(pow(r,2)*pi)
            l_ofs.append(r)    
        else:
            r = 20
            l_ofs.append(r)
            rs.append(pow(r,2)*pi)

        for c in children:
            x1,y1 = node_xs[c], node_ys[c]
            p = patches.ConnectionPatch([x0,y0],[x1,y1],'data','data',
                                        arrowstyle = 'wedge,tail_width=3',
                                        shrinkA = 25,
                                        shrinkB = 25,
                                        edgecolor = 'black',
                                        facecolor = 'white',
                                        alpha = .5)
            cxns.append(p) 
        
    f = plt.figure(0)
    f.clear()
    ax = f.add_subplot('111')
    for x in cxns:
        ax.add_patch(x)
        
    ax.scatter(xs,ys,rs,'white')
    for i in range(len(xs)):
        ax.annotate(ls[i],[xs[i],ys[i]],
                    xytext = [l_ofs[i],l_ofs[i]],
                    textcoords = 'offset points')
        
