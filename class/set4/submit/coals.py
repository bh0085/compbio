from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import logical_and as land
from numpy import logical_or as lor

def run0(n = 100):
    n10 =n/2    

    nruns = 1000
    nrep = 1000
    ys = zeros((nrep, nruns)) -1
    lifetimes = zeros(nrep)

    for i in range(nrep):
        n1 = n10
        for j in range(nruns):
            ys[i,j]= n1
            frac = float(n1) / n
            r = random.uniform(0,1,n)
            n1 = len(nonzero(less(r,frac))[0])

            if n1 == n or n1 == 0:
                ys[i,j:] = n1
                lifetimes[i] = j
                break
    
    
    f =plt.figure(0)
    f.clear()
    ax = f.add_subplot(211)
    for i in range(nrep):
        ax.plot(ys[i])
        
    ax2 = f.add_subplot(212)
    ax2.hist(log(lifetimes))
        

def run_meta():
    ntrials = 1000
    nk = 3
    means = zeros((ntrials,nk))
    stds = zeros((ntrials,nk))

    
    f = plt.figure(0)
    f.clear()

    for i in range(ntrials):
        m,s = run(n = 500, gens = 3500)
        means[i,:] = m
        stds[i,:] = s
        print i

    ax1 =  f.add_subplot('211')
    ax2 = f.add_subplot('212')
    
    print 'means:'
    for m in means.T:
        ax1.plot(m)
        print mean(m)
    print 'stds:'
    for s in stds.T:
        ax2.plot(s)
        print mean(m)
    return means, stds

       

def run_sex_meta():
    ntrials = 1000
    nk = 3
    means = zeros((ntrials,nk))
    stds = zeros((ntrials,nk))

    for i in range(ntrials):
        m,s = run_sex()
        means[i,:] = m
        stds[i,:] = s
        print i
    
    print 'means:'
    for m in means.T:
        print mean(m)
    return means


def run1(n = 100, gens= 3000, plot_pop = 0):

    pop = [[i] for  i in range(n)]
    for t in range(gens):
        ids = array(floor(random.uniform(0,n,n)),int)
        newpop = []

        for i in range(n):
            row = list(pop[ids[i]])
            row.append(i)
            newpop.append(row)
        pop = newpop
        
    pop = array(pop,float)
            
    if plot_pop:
        f = plt.figure(0)
        ax = f.add_subplot(211)
        ax.imshow(pop/n,interpolation = 'nearest',aspect = 'auto')

    


    plot_total = 0
    if plot_total:
    #step backward in time to count coal:
        unq = zeros(n,int) + 1
        c_t = zeros(gens)
        for t in range(gens)[::-1]:
            h = zeros(n)
            col = pop[:,t]
            for elt in col[nonzero(unq)[0]]:
                h[elt]+=1

            collisions = nonzero(h > 1)[0]
        
            for c in collisions:
                marks = nonzero(land(unq,equal(col,c)))[0][1:]
                c_t[t] += len(marks)
                unq[marks] = 0
    

        ax2 = f.add_subplot(212)
        ax2.plot(c_t[::-1][:100])
            
    plot_k = 1
    kvals = [2,3,4]
    ntrials = 100
    ctimes =zeros((len(kvals),ntrials),int)
    for j in range(len(kvals)):
        k = kvals[j]
        for i in range(ntrials):
            entries = array(random.uniform(0,n,k),int)
            subarr = pop[entries,:]
            ctime = -1
            for t in range( gens )[::-1]:
                if len(nonzero(not_equal(subarr[:,t], subarr[0,t]))[0]) == 0:
                    ctime = gens -t
                    break
            ctimes[j,i] = ctime

    mean_ctimes = []
    std_ctimes = []
    for k in range(len(kvals)):
        mean_ctimes.append(  np.mean(ctimes[k,:]))
        std_ctimes.append(std(ctimes[k,:]))
    return  mean_ctimes, std_ctimes


    

    
def run_sex(n = 500, gens =3500, F = 100):
    M = n-F
    mdist = random.randint(0,M,(n,gens)) + F
    fdist = random.randint(0,F,(n,gens))
    which = random.randint(0,2,(n,gens))
    pop = fdist
    nz = nonzero(which)
    pop[nz] = mdist[nz]
    
    
    plot_k = 1
    kvals = [2,3,4]
    ntrials = 1
    ctimes =zeros((len(kvals),ntrials),int)
    for j in range(len(kvals)):
        k = kvals[j]
        for i in range(ntrials):
            entries = array(random.uniform(0,n,k),int)
            for ct in range(gens):
                entries = pop[entries,ct]
                if var(entries)  == 0: break
                if ct == gens -1:
                    print 'SHIT! NEVER COALESCED'
            ctimes[j,i] = ct

    mean_ctimes = []
    std_ctimes = []
    for k in range(len(kvals)):
        mean_ctimes.append(  np.mean(ctimes[k,:]))
        std_ctimes.append(std(ctimes[k,:]))
    return  mean_ctimes, std_ctimes
    


                    
def run(n = 500, gens = 5000):
    pop = array(floor(random.uniform(0,n,(n,gens))),int)
       
    plot_k = 1
    kvals = [2,3,4]
    ntrials = 1
    ctimes =zeros((len(kvals),ntrials),int)
    for j in range(len(kvals)):
        k = kvals[j]
        for i in range(ntrials):
            entries = array(random.uniform(0,n,k),int)
            for ct in range(gens):
                entries = pop[entries,ct]
                if var(entries)  == 0: break
                if ct == gens -1:
                    print 'SHIT! NEVER COALESCED'
            ctimes[j,i] = ct

    mean_ctimes = []
    std_ctimes = []
    for k in range(len(kvals)):
        mean_ctimes.append(  np.mean(ctimes[k,:]))
        std_ctimes.append(std(ctimes[k,:]))
    return  mean_ctimes, std_ctimes
    
