import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import compbio.utils.colors as mycolors

def test_kmeans(seed = 0):
    k =2
    n0 = 200
    n1 = 400
    m0 = array([[20.,0],
                [0,30.]])
    v0 = array([[10,20],
                [20,10]])
    
    random.seed(seed)
    ndim = 4
    
    freqs = zeros((2,ndim),float) + .1
    freqs[0,0] = .9
    freqs[0,2] = .5
    freqs[1,1] = .9
    freqs[1,3] = .5
    data = zeros((n0+n1,ndim))

    
    
    for idx in range(ndim):
        data[:n0,idx] = random.binomial(1,freqs[0,idx],n0)
        data[n0:,idx] = random.binomial(1,freqs[1,idx],n1)
        
    print data[0,:]
    #start kmeans with default values
    test_withdata(k,data)

def test_withdata(k,data):
    km = kmeans(k, data)
    km.fig=10
    km.draw()
    while 1:
        km.run_n(10)
        km.fig=11
        km.drawspec()
        km.fig=12
        km.draw()
        inp = raw_input('done? (d)')
        if inp =='d':
            break

class kmeans():



    def __init__(self,k, data, 
                 initialization = None,
                 seed = 0,
                 prior = 'prop',
                 max_runs = 100):
        
        print 'WHAT THE HELL'
        self.fig = 0
        self.max_runs = max_runs

    
        self.data = data

        sd = shape(data)
        self.nd = sd[0]
        self.labels = zeros(self.nd,int)

        self.dim = sd[1]
        self.nm = k
        self.max_runs = max_runs

        #RANDOM initialization
        if not initialization:
            #data_mins = np.min(data,0)
            #data_maxs = np.max(data,0)

            random.seed(seed)
            self.means = random.uniform(0,1,(self.nm,self.dim))

        else:
            raise Exception('non random initialization not yet implemented')

        #get a prior on each dimension from the data (or flat)
        if prior == 'prop':
            normalizer = np.array(np.sum(self.data,0),float)
            normalizer[nonzero(equal(normalizer, 0))[0]] = 1
            self.priors = normalizer
        elif prior =='flat':
            self.priors = 1+np.zeros(self.dim) 
        self._compute_labels()

    def run_n(self, n):
        for i in range(n):
            self._iterate()
            

    def run_till_convergence(self,
                             tol=.05):
        scale = np.mean(np.sum(np.pow(data,2),1),0)
        runs = 0 
        while 1:
            runs += 1
            means_0 = array(self.means)
            self._iterate()
            divergence = np.mean(np.sum(np.pow(means_0 - means,2),1),0)
            if divergence / scale < tol:
                failed = False
                break
            if runs > self.max_runs:
                failed = True
                break
            
        return failed
    def draw2d(self):
        f = plt.figure(self.fig)
        f.clear()
        ax = f.add_axes([.05,.05,.9,.9])
        data_x = 0
        data_y = 1
        
        ct = mycolors.getct(self.nm)
        xs, ys, rs, cs = [[] for i in range(4)]
        for i in range(self.nd):
            xs.append(self.data[i][data_x])
            ys.append(self.data[i][data_y])
            rs.append(25)
            cs.append(ct[self.labels[i]])
            
        for i in range(self.nm):
            xs.append(self.means[i][data_x])
            ys.append(self.means[i][data_y])
            rs.append(100)
            cs.append([0,0,0])
        ax.scatter(xs,ys,rs,cs)
    
    def draw(self):
        f = plt.figure(self.fig)
        f.clear()
        ax = f.add_axes([.05,.05,.9,.9])
        for m in self.means:
            ax.plot(m)

    def drawspec(self):
        f = plt.figure(self.fig)
        f.clear()
        img = self.means
        mean_idxs = arange(self.nm)
        weights = array(np.sum(equal(mean_idxs[:,newaxis],self.labels[newaxis,:]),1),float)
        
        weights /= (sum(self.means * self.priors,1))

        ax = f.add_axes([.05,.05,.9,.9])
        activity = mycolors.imgcolor(img,BW = True)
        ax.imshow(activity,aspect = 'auto')

        activity = mycolors.imgcolor(img*weights[:,newaxis],alpha = True)
        ax.imshow(activity,aspect = 'auto',interpolation = 'nearest')
        
    def _iterate(self):
        #compute labeling of data points given means
        self._compute_labels()
        #given labeling, compute new means
        self._compute_means()
            
    
    def _compute_labels(self):

        #compute labels from the max of the log likelihood of each data point
        #under a gaussian with width_d = prior[d]

        self.labels = np.argmax(np.sum( self.data[:,newaxis,:] * log(self.means/self.priors),2),1)

        #for i in range(self.nd):
        #    self.labels[i] = argmax([ np.sum(self.data[i,:] * log (self.means[j,:]
        #                                                           / self.priors )
        #                                     )
        #                                      
        #                            for j in range(self.nm)
        #                            ])
    def _compute_means(self):
        
        #recompute means from the membership of each cluster
        for i in range(self.nm):
            elts = nonzero(equal(self.labels,i))[0]
            #for nwo, we just don't bother updating when a mean 
            #has no elements.
            alpha = .2
            if len(elts) != 0:
                self.means[i] = (1-alpha)*self.means[i] + alpha*mean(self.data[elts,:],0)
            

