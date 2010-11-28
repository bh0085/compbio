from numpy import *
import numpy as np
import tree

m = [[0,3,11,10,12],
     [3,0,12,11,13],
     [11,12,0,9,11],
     [10,11,9,0,8],
     [12,13,11,8,0]]

class tree_alg():
    def __init__(self,pairwise):
        self.pairwise = pairwise
        self.tree = None
    def compute(self, max_itr = -1):
        dists = array(self.pairwise)
        
        n = len(dists)
        self.T = range(n)
        self.L =self.T
        parents = [-1 for i in range(n)]
        children = [[] for i in range(n)]
        weights = [[] for i in range(n)]
        
        itr = 0 
        while 1:
            itr+= 1
            if itr == max_itr: break
            if len(self.L) == 2:
                parents.append(-1)
                new_idx = len(self.T)
                for l in self.L:
                    parents[l] = new_idx
                children.append(self.L)
                weights.append([dists[0,1]/2,dists[0,1]/2])
                dists = []
                break
                

            #compute custom distances of nj...
            self.computeD(dists)
            #get a minimum pair according to our alg.
            lilidx = self.min_pair(dists)
            #compute distances to our new node.
            lildists = self.new_dist(dists,lilidx)


            nd = len(dists)
            keep = range(nd)
            for i in lilidx: keep.remove(i)

            new_idx = len(self.T)
            self.T.append(new_idx)

            #compute edge weights
            print 'computing edge weights before adding in the new node... check to see if this gives right answer'
            weights.append(self.edgeweights(dists,lilidx))


            parents.append(-1)
            children.append(map(lambda x: self.L[x],lilidx))
            for i in lilidx :parents[self.L[i]] = new_idx

            self.L = list(array(self.L)[keep])
            self.L.append(new_idx)
            
            if nd <=2: break
            lildists = lildists[keep]
            dists_new = zeros((nd -1,nd-1))

            dists_new[:nd-2,:nd-2] = dists[keep][:,keep]
            dists_new[nd-2,:nd-2]= lildists
            dists_new[:nd-2,nd-2]= lildists
            dists = dists_new

        labels = ['' for i in range(len(parents))]
        labels[0:5] = ['a','b','c','d','e']
        self.tree = tree.Tree2(parents,children,weights = weights, labels = labels)
        return dists

class upgma(tree_alg):

    def computeD(self,dists):
        self.r = zeros(len(dists))
        self.D = dists

    def min_pair(self,dists):
                    #find the closest pair
        lil = np.min(dists[nonzero(greater(dists,0))])
        lilidx = nonzero(equal(dists,lil))
        lilidx = map(lambda x:x[0],lilidx)
        return lilidx

    def new_dist(self,dists,pair):
        davg = mean(dists[pair,:],0)
        return davg
    def edgeweights(self,dists,lilidx):
        return [1,1]

class nj(tree_alg):
    def computeD(self,dists):
        n = len(dists)
        r = zeros(n,float)
        D = zeros((n,n),float)
        L = self.L
        for i in range(n):
            r[i]= 1.0/ (len(L) -2) * np.sum(dists[i])
        for i in range(n):
            for j in range(n):
                if i ==j: continue
                D[i,j] = dists[i,j] -(r[i] + r[j])
        
        self.r = r
        self.D = D

    def min_pair(self,dists):
                    #find the closest pair
        lil = np.min(self.D)
        lilidx = nonzero(equal(self.D,lil))
        lilidx = map(lambda x:x[0],lilidx)
        return lilidx

    def new_dist(self,dists,pair):
        l = len(dists)
        dnew = zeros(l)
        for i in range(l):
            dnew[i] = .5 * (dists[pair[0],i] +
                            dists[pair[1], i] -
                            dists[pair[0],pair[1]])

        return dnew
    def edgeweights(self,dists,pair):
        return [.5 * (dists[pair[0],pair[1]] + self.r[pair[0]] - self.r[pair[1]]),
                 .5 * (dists[pair[0],pair[1]] + self.r[pair[1]] - self.r[pair[0]])]
        


        
def get_samples(rnd = False):
    if rnd:
        dim = 2
        n = 10
        vecs = random.uniform(0,5,[n,dim])
        diffs = sum(pow(vecs[:,:],2),1)[:,newaxis] + \
            sum(pow(vecs[:,:],2),1)[newaxis,:] + \
            - sum(2*vecs[:,newaxis,:]*vecs[newaxis,:,:],2)
        return diffs
    else:
        return m



def run_test_upgma(max_itr = -1):
    s = get_samples()
    u = upgma(s)
    dists = u.compute(max_itr)
    return u,dists


def run_test_nj(max_itr= -1):
    s = get_samples()
    u = nj(s)
    dists = u.compute(max_itr)
    return u, dists





