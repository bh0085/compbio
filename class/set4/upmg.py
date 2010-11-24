from numpy import *
import numpy as np

class UPMG():
    def __init__(self,pairwise):
        self.pairwise = pairwise
    def compute(self):
        dists = array(self.pairwise)
        
        n = len(dists)
        nodes = range(n)
        node_idxs = nodes

        parents = [-1 for i in range(n)]
        children = [[] for i in range(n)]

        while 1:
            #find the closest pair
            lil = np.min(dists[nonzero(greater(dists,0))])
            lilidx = nonzero(equal(dists,lil))
            lilidx = map(lambda x:x[0],lilidx)
            
            #average together the distance vectors.
            davg = mean(dists[lilidx,:],0)
            nd = len(dists)
            keep = range(nd)
            for i in lilidx: keep.remove(i)

            new_idx = len(nodes)
            nodes.append(new_idx)

            parents.append(-1)
            children.append(lilidx)
            for i in lilidx :parents[node_idxs[i]] = new_idx

            node_idxs = list(array(node_idxs)[keep])
            node_idxs.append(new_idx)

 
            
            if nd <=2: break
            davg = davg[keep]
            dists_new = zeros((nd -1,nd-1))

            dists_new[:nd-2,:nd-2] = dists[keep][:,keep]
            dists_new[nd-2,:nd-2]= davg
            dists_new[:nd-2,nd-2]= davg
            dists = dists_new

        print parents

def get_samples():
    dim = 2
    n = 4
    vecs = random.uniform(0,5,[n,dim])
    diffs = sum(pow(vecs[:,:],2),1)[:,newaxis] + \
        sum(pow(vecs[:,:],2),1)[newaxis,:] + \
        - sum(2*vecs[:,newaxis,:]*vecs[newaxis,:,:],2)
    return diffs



def run():
    s = get_samples()
    u = UPMG(s)
    u.compute()

