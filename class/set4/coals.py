from numpy import *
import numpy as np
import matplotlib.pyplot as plt

def run(n = 100):
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
        
