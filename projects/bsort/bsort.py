from numpy import *
import numpy as np
import compbio.utils.plots as myplots
import matplotlib.pyplot as plt
import scipy.signal as ss

def sample_array(nx = 100, ny = 100, t = 'default', perc_high = .02):

    
    if t == 'default':
        dims = (nx,ny)
    elif t== 'RGB':
        dims = (nx,ny,3)

    uval = random.uniform(0,.25,dims)
    ntot = nx * ny
    nout = ntot*perc_high
    if t=='RGB': nout*= 3
    if t=='RGB':
        uval[random.random_integers(0,nx -1,nout),
             random.random_integers(0,ny -1,nout),
             random.random_integers(0,2,nout)]=random.uniform(.25,1,nout)
    else:
        uval[random.random_integers(0,nx -1,nout),
             random.random_integers(0,ny -1,nout)]=random.uniform(.25,1,nout)
        
    return uval
def run0(meth = 'moment', 
         nx = 20, 
         ny = 20, 
         fig = 1, 
         seed = 0,
         nfl = 200,
         itr = 2,
         arr = None,
         transpose = False):
    '''Generates a random array and clusters red values along the diagonal.'''
    random.seed(seed)
    if arr == None:
        s = sample_array(nx = nx, ny = ny,t = 'default',perc_high = .2)
    else:
        nx, ny  = shape(arr)
        s = arr

    qbad = 0

    all_scores = []
    all_dlts = []
    for i in range(itr):
        
        temp = (1 -  float(i)/itr)
        dims = shape(s)
        score = makescore(temp, dims)
        all_scores.append(sorted(reshape(score,-1)))
        #if qbad:
            #score = score.T
        
        if meth == 'moment':
            s = get_moment_sort(s)
        if meth == 'maxfirst':
            s = get_maxfirst(s)
        elif meth == 'rnd':
            this_dlt = []
            s = rnd_shuffles(s,
                             nfl = nfl * (float(shape(s)[0])/max(shape(s))), 
                             score = score,
                             temp = temp,
                             deltas = this_dlt)
            all_dlts.append(this_dlt)

        if transpose:
            qbad = 1-qbad
            s = s.T
            #score = score.T

        #if i > 10:

    #s = ss.medfilt2d(s,5)

    f2 = plt.figure(fig)
    f2.clear()
    ax = f2.add_subplot(111)
    ax.imshow(s[:,:,newaxis]*[1,0,0], interpolation = 'nearest')
    
    f3 = plt.figure(fig+1)
    f3.clear()
    ax2=f3.add_subplot(211)
    for s in all_scores:
        ax2.plot(s)
    
    ax3 = f3.add_subplot(212)
    for s in all_dlts:
        deltas = array(map(lambda x: x[0],s))
        acc = array(map(lambda x: x[1],s))
        ax3.plot(sorted(deltas))
    

def makescore(temp,dims):
    score = zeros(dims) + 1
    
    lonax = greater(dims[1],dims[0])
    londim= dims[lonax]
    sdim = dims[1-lonax]
    ldiag = arange(0,londim,1,int)
    sdiag = array(floor(ldiag * (float(sdim)/londim)),int)
    dvals = [ldiag, sdiag]
    if lonax:dvals = dvals[::-1]
    score[dvals] = 20.*temp ** .5

    if lonax: score = score.T

    #for i in range(len(score)):
    #    score[i][i] = 20.* temp **.5

    if temp > .1:
        g = ss.gaussian((londim/2)*temp**2,(londim/2)*temp**2)[:,newaxis]*\
            ss.gaussian((sdim/2)*temp**2,(sdim/2)*temp**2)[newaxis,:]

        g/= sum(g)
        score = ss.convolve2d(score,g,'same','wrap')
        
    for i in range(1,len(score)):
        score[i,arange(sdiag[i])] *= .25


        
    if lonax: score = score.T

    return score

def get_moment_sort(arr):
    nx, ny= shape(arr)
    yvals = array((reshape(arange(nx *ny), (nx, ny)) / nx),float)
    moments = np.sum(yvals[:,:] * arr,0)/np.sum(arr,0)
    srt_red = argsort(moments[:])[::-1]
    return arr[:,srt_red]

def get_maxfirst(arr):
 
    dims = shape(working_arr)
    wal = reshape(working_arr,-1)

    sdim = dims[0]
    ddim = dims[1]

    biggest = argsort(wal)[::-1]
    diag_used = zeros(ddim)
    shuf_used = zeros(sdim)
    shuf_tgs = []
    for b in biggest:

        if len(shuf_tgs)> min(dims): break

        bs,bd = divmod(b,dims[1])
        if not(diag_used[bd] or shuf_used[bs]):
            diag_used[bd] =1
            shuf_used[bs] =1
            shuf_tgs.append([bs,bd])

    out = zeros(dims)
    for i in range(len(shuf_tgs)):
        out[shuf_tgs[i][1]] = working_arr[shuf_tgs[i][0]]
        
    return out

def rnd_shuffles(arr, nfl = 100,qbad = 0,score = None,temp = None, deltas = None):
    '''
    Use random shuffles to maximize the scoring function.
    '''


    print
    print shape(score), shape(arr)
    print nfl
    print

    record_deltas = False
    if deltas != None:
        record_deltas = True
        
    dims = shape(arr)
    t_norm = std(score)
    
    for i in range(nfl):
        pair = random.random_integers(0,dims[0]-1,2)
        rows = array(arr[pair, :])
        s0 = np.sum(score[pair,:] * rows)
        s1 = np.sum(score[pair[::-1],:] *rows)
        diff = s1 -s0
        if s1 > s0:
            arr[pair[1]] = rows[0]
            arr[pair[0]] = rows[1]
            acc = True
        else:
            delta =(s0 -s1)/(t_norm)
            p = exp(-delta / (temp/4))
            if random.uniform(0,1) < p:
                arr[pair[1]] = rows[0]
                arr[pair[0]] = rows[1]
                acc = True
            else:
                acc = False
                
        if record_deltas:
            deltas.append([diff,acc])
            
    
    return arr
                    
