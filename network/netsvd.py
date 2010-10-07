import pickle
from numpy import *
import matplotlib.pyplot as plt
def compute(name, number):
    #possible names: bRN, kRN, fRN, mRN
    trgs, tfs = nu.parse_net(name)
    print 'N Targets: ' + str(len(trgs.keys()))
    print 'N TFs: ' + str(len(tfs.keys()))
    
    k0 = trgs.keys()
    k1 = tfs.keys()

    kmap0 = {}
    kmap1 = {}
    for i in range(len(k0)):
        kmap0[k0[i]] = i
    for i in range(len(k1)):
        kmap1[k1[i]] = i

    #Now, for a first shot, lets construct an n*m matrix
    #with possible interactions for each gene appearing in 
    #the network        
    #prune the network for testing purposes.
    m = len(kmap0.keys())

    n = len(kmap1.keys())
    N = np.zeros([m,n])
    for k_trg,v_trg in trgs.items():
        i = kmap0[k_trg]
        for k_tf in v_trg['tfs']:            
            j = kmap1[k_tf]
            N[i][j] = 1.0

    import pickle
    svdname = name + '.svd'
    comp_svd = True
    if comp_svd:
        U, S, Vh = lin.svd(N)
        V = Vh.T
    else:
        U,S,V = pickle.load(open(svdname))

    #from mlabwrap import mlab
    #U2, S2, V2 = mlab.svd(N,nout = 3)
        
    V = mat(V)

    uvecs = U[:,:len(S)]
    vvecs = V
    dosave = ((kmap0,kmap1),(uvecs, S, vvecs))
    pickle.dump(dosave,open('temp/'+name+'_svd.pickle','w'))

def view_net(name = 'mRN'):
    if type(name) == type(''):
        name = [name]
    else:
        if 'all' in name:
            name =  ['mRN','fRN','kRN','bRN']
    for i in range(len(name)):
        net_name = name[i]
        keys, svd = pickle.load(open('temp/'+net_name+'_svd.pickle'))
        U, S, V = svd
        view_svd(U,S,V, fig = i)

def view_svd(U,S,V, fig = 0 ):    
    m = shape(U)[0]
    n = shape(V)[0]

    nt = shape(U)[0]
    yvals = pow(reshape(mod(arange(m*n) , n),(m,n)),2)
    tots = sum( abs(U) *yvals,1)
    srt = argsort(tots)
    Usrt = U[srt,:]
        
    fig = plt.figure(fig)
    fig.clear()
    ax = fig.add_axes([0,0,1,1],frameon = False)
    ax.imshow(abs(U[srt[0:m/2],:]),aspect = 'auto',interpolation = 'nearest')

