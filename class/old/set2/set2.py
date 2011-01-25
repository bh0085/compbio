from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import random
import itertools as it
import pickle

state_seed = 0

def random_rna(n,gc = .5, seed = 1):
    random.seed(seed)
    seq = ''
    for i in range(n):
        r = random.uniform(0,1)
        if r < gc:
            if random.uniform(0,1) > .5:
                seq+='C'
            else:
                seq+='G'
        else:
            if random.uniform(0,1) > .5:
                seq+='A'
            else:
                seq+='U'    
    return seq

def run(part = '1a',rseed = 1):
    sequence = 'MLAN'
    if part == '1a':
        T,E,skey,ekey = init_protein()
        draw_hmm(T,E,skey,ekey)
    if part == '1b':
        T,E,skey,ekey = init_protein()
        M, ptrs = viterby_path(T,E,skey,ekey,sequence=sequence,start_state='Other')
        path = viterby_traceback(M,ptrs)
        draw_path(path,sequence,M,T,E,skey,ekey)
    if part =='1d':
        T,E,skey,ekey = init_protein()
        M,N = viterby_sum(T,E,skey,ekey,sequence=sequence,start_state='Other')
        draw_sum(M,N,T,E,skey,ekey)
    if part =='2a':
        M1,M2 = init_p2()
        seqs1, seqs2 = [], []
        for i in range(200):
            seqs1.append( p2_randomseq(M1,6,i+9000))
            seqs2.append( p2_randomseq(M2,6,i+2000))
            
        
        p1m1 , p1m2 = [],[]
        for s in seqs1:
            p1m1.append(p2_prob( M1,s))
            p1m2.append( p2_prob( M2,s))
        
    
        p2m1 ,p2m2 = [],[]
        for s in seqs2:
            p2m1.append(p2_prob(M1,s))
            p2m2.append(p2_prob(M2,s))
        p2m1 = array(p2m1)
        p2m2 = array(p2m2)
        p1m1 = array(p1m1)
        p1m2 = array(p1m2)
        
        p2_drawprobs( (p1m1,p1m2),(p2m1,p2m2) )
    if part == '4':
        r = random_rna(50)
        r = 'CGCUCAUCCAGAGGGGCAGAGGGAUACGGCCCGAUGAAGCCCCGGCAACCCUCCAGUCGGUUCUUGUCACACGGACGUGGCGAGGCUCCCGGCUAGGGAAGGUGCCAAAUCCGUCUCACGGCGAGAUGCGUCGUGAGGAAGAUGAGGA'
        #S. coelicolor
        bac ='CUCUUAUCAAGAGCAGGCAGAGGGACUUGGCCCGAUGAUGCCCAGCAACCGACCGUAAUUCCAUCGUGAGAUGGGGCGCAAUCCUUGCGCCGGAGAAUUCCUCCAUAAGGCACGGUGCUAAUUCCAGCAGAAAGCUUGGCUUUCUGGCAGAUAAGAG'
        r =bac
        M,TB = nussinov(r)
        t = nussinov_traceback(r,M,TB)
        draw_rna(r,t)
    if part == '4a':
        scores = []
        for i  in range(200):
            if not mod(i,10): print i
            seq = random_rna(100,gc = .5, seed=rseed)
            M,TB = nussinov(seq)
            score = (M[0,-1])
            scores.append(score)
        scores = array(scores)
        print mean(scores)
    if part == '4b':
        get_scores = False
        if get_scores:
            all_scores, all_is = [], []
            for i in np.power(arange(4,17),2):
                print 'computing, i = ' + str(i)
                scores= []
                for j in range(20):
                    seq = random_rna(i,gc = .5, seed=rseed)
                    M,TB = nussinov(seq)
                    score = (M[0,-1])
                    scores.append(score)
                mscore = mean(array(scores))
                all_is.append(i)
                all_scores.append(mscore)
            scores = {'scores':all_scores,
                      'l':all_is}
            pickle.dump(scores,
                        open('4b_scores.pickle','w'))
        else:
            scores = pickle.load(open('4b_scores.pickle'))
        p4_plot(scores)
    if part == '4c':
        get_scores = False
        if get_scores:
            all_scores, all_is = [], []
            for i in linspace(.01,.99,100):
                print 'computing, i = ' + str(i)
                scores= []
                for j in range(20):
                    seq = random_rna(50,gc = i, seed=rseed)
                    M,TB = nussinov(seq)
                    score = (M[0,-1])
                    scores.append(score)
                mscore = mean(array(scores))
                all_is.append(i)
                all_scores.append(mscore)
            scores = {'scores':all_scores,
                      'l':all_is}
            pickle.dump(scores,
                        open('4c_scores.pickle','w'))
        else:
            scores = pickle.load(open('4c_scores.pickle'))
        p4_plot(scores)
            
                
                


def getct(n, seed = 0, forceRB = True):
    random.seed(seed)
    ct = []
    for i in range(n):
        if forceRB:
            if i == 0:
                ct.append([1,0,0])
                continue
            elif i == 1:
                ct.append([0,0,1])
                continue
        ct.append([random.uniform(0,1),
                   random.uniform(0,1),
                   random.uniform(0,1)])
    return ct



def viterby_path(T,E,skey,ekey,sequence = 'MLAN',start_state = 'Other'):
    l = len(sequence)
    nt = shape(T)[0]
    M = zeros((nt,l))
    ptrs = zeros((nt,l))
    for j in range(l):
        for i in range(nt):
            if j == 0:
                if i == skey[start_state]:
                    tprob = 1.0
                else:
                    tprob = 0.0
                max_idx =skey[start_state]
            else:
                max_idx =argmax( M[:,j-1]*T[:,i] )
                tprob=( M[:,j-1]*T[:,i] )[max_idx]
            e = ekey[sequence[j]]
            eprob = E[e,i]
            M[i,j] = eprob*tprob
            ptrs[i,j] = max_idx

    return (M,ptrs)
def viterby_sum(T,E,skey,ekey,sequence = 'MLAN',start_state = 'Other'):
    l = len(sequence)
    nt = shape(T)[0]
    M = zeros((nt,l))
    N = zeros((nt,l))

    ptrs = zeros((nt,l))
    for j in range(l):
        for i in range(nt):
            if j == 0:                
                if i == skey[start_state]:
                    tprob = 1.0
                    npaths = 1
                else:
                    tprob = 0.0
                    npaths = 0
                max_idx =0
                
            else:
                #max_idx =argmax( M[:,j-1]*T[:,i] )
                tprob=sum( M[:,j-1]*T[:,i] )
                path_inds = nonzero(T[:,i])
                npaths = sum( N[path_inds,j-1])
            e = ekey[sequence[j]]
            eprob = E[e,i]
            if eprob > 0:
                N[i,j] = npaths
            M[i,j] = eprob*tprob
            ptrs[i,j] = max_idx
    return (M,N)

def viterby_traceback(M,ptrs):
    print ptrs
    print M
    ns,l = shape(M)
    max_path = []
    
    max_end = argmax(M[:,-1])
    ptr = ptrs[max_end,l-1]
    max_path = [max_end]
    cur = l - 2
    while cur >= 0:
        next = ptr        
        max_path.append(ptr)
        ptr = ptrs[next,cur]        
        cur-=1

    return max_path[::-1]

def init_protein():
    T = array([[.7,.1,.2,.0],
               [.2,.6,.0,.2],
               [.3,.3,.1,.3],
               [.3,.3,.0,.4]])

    E = array([[.35,.10,.00,.05],
               [.30,.05,.15,.15],
               [.15,.30,.10,.20],
               [.10,.40,.10,.15],
               [.05,.00,.35,.20],
               [.05,.15,.30,.25]])
         
    skey={'Alpha Helix':0,
          'Beta Sheet':1,
          'Turn':2,
          'Other':3}
    ekey={'M':0,
          'L':1,
          'N':2,
          'E':3,
          'A':4,
          'G':5}
    return(T,E,skey,ekey)

def draw_hmm(T, E, skey, ekey):
    syms = sym_map()
    items = list(skey.itervalues())
    smap = [  syms[skey.keys()[items.index(i)]] for i in range(max(items)+1)  ]
    items = list(ekey.itervalues())
    emap = [  syms[ekey.keys()[items.index(i)]] for i in range(max(items)+1)  ]


    st = T.shape
    se = E.shape
    nt = st[0]
    ne = se[0]
    
    f = plt.figure(0)
    f.clear()
    max_x = 10
    max_y = 10
    ax = f.add_axes([0,0,1,1],xlim = [-2,max_x+2],ylim = [-2,max_y+2])
    ax.autoscale(False)

    scolors = getct(nt,seed = state_seed)
    scolor = 'blue'
    ecolor = 'red'
    falpha = .7
    
    node_r = 6000
    sx, sy, ex, ey = [], [], [], []


    
    for i in range(nt):
        x = (i+1) * float(max_x)/(nt +1)
        y = max_y
        r = node_r
        sx.append(x)
        sy.append(y)
        
        ax.scatter(x,y,r*1.5, c = 'white', edgecolor = 'none',zorder = i + ne - .1)
        ax.scatter(x,y,r, c ='none',
                   edgecolor = scolors[i],linewidth = 15,
                   alpha = falpha, zorder = i+ne)

        key = smap[i]
        tsize = 20
        xofs = 10*( i - (ne-2)/2)
        ax.annotate(key,xy = (x,y), xycoords = 'data',
                    size=tsize,
                    va = 'center',
                    ha = 'center',
                    zorder = nt+ne,
                    bbox=dict(boxstyle="round", fc="1",alpha = .6,ec='none'))
        
    for i in range(ne):
        x = (i+1) * float(max_x) /(ne+1)
        y = 0
        r = node_r/2
        ex.append(x)
        ey.append(y)
        ax.scatter(x,y,r*1.5,c = 'white', edgecolor = '0',zorder = -1)
        #ax.scatter(x,y,r,c = ecolor, alpha = .1)

        key =emap[i]
        tsize = 20
        xofs = 10*(i - ne/2)
        ax.annotate(key,xy = (x,y),xycoords = 'data',
                    va = 'center',
                    ha = 'center',
                    size=tsize,
                    zorder = nt+ne,
                    bbox=dict(boxstyle="round", fc="1",alpha = .6,ec = 'none'))
        
    for i in range(nt):
        max_prob = max(T[i,:])
        for j in range(nt):
            if i == j: continue
            wt = T[i,j]
            
            p = T[i,j]
            aalpha = p
            if p >= max_prob*.75: aalpha = 1
            
            hwidth = 1

            cur_z = ne +  nt + 1

            a_mag = 90
            aofs = 90 - a_mag/2
            cstr = 'arc3,rad = .4'

            ax.annotate('', xy=(sx[i]+.1,sy[i]+.1),  xycoords='data',
                        xytext=(sx[j],sy[j]), textcoords='data',
                        size=20,
                        bbox=dict(boxstyle="round", fc="0.8"),
                        arrowprops=dict(arrowstyle="fancy",
                                        fc="0", ec="none",
                                        #patchB=el,
                                        alpha = aalpha,
                                        connectionstyle=cstr,
                                        relpos = [.3,-.3]
                                        ,shrinkA = 50,shrinkB=70),
                        zorder = cur_z
                        )
        
    for i in range(st[0]):
        for j in range(se[0]):
            x,y = sx[i],sy[i]
            xf, yf = ex[j],ey[j]
            dx, dy = xf - x , yf - y

            max_prob = max(E[j,:])
            p = E[j,i]
            aalpha = p
            if p >= max_prob*.75: aalpha = 1

            hwidth = 1

            astr = ''
            if wt > .3:
                astr = str(wt)
            
            a_mag = 30
            aofs = -90 - a_mag/2
            #cstr = "angle3,angleA="+str(float(j+1)/(ne+1)*a_mag +aofs)+",angleB="+str(float(i+1)/(nt+1)*0 + 0)
            cstr = "arc3,rad="+str(-(i - nt/2) *.2)
                
            ax.annotate('', xy=(xf,yf),  xycoords='data',
                        xytext=(x,y), textcoords='data',
                        size=20,
                        bbox=dict(boxstyle="round", fc="0.8"),
                        arrowprops=dict(arrowstyle="fancy",
                                        fc="0", ec="none",
                                        #patchB=el,
                                        connectionstyle=cstr,
                                        alpha = aalpha,zorder = -5,
                                        shrinkA=20,shrinkB = 20),zorder = -5
                        )

    fakeax = f.add_axes([0,0,1,1],visible = False)
    legs = (fakeax.plot(1,1,markeredgecolor = '0',markerfacecolor = [1,0,0],
                     visible = True,linestyle = 'none',
                     marker = '.',markersize = 20,label = 'Emission'),
            fakeax.plot(1,1,markeredgecolor ='0',markerfacecolor = [0,0,1],
                     visible = True,linestyle = 'none',
                     marker = '.',markersize = 20,label = 'State'),
            fakeax.plot(1,1,color = '0',alpha=1,
                     visible = True,linestyle = '-',linewidth = 3,
                     label = 'Probable \nEmission/Transition'),
            fakeax.plot(1,1,color = '0',alpha=.3,
                     visible = True,linestyle = '-',linewidth = 3,
                        label = 'Improbable \nEmission/Transition \n(Opacity gives likelihood)'),
            fakeax.plot(1,1,color = '0',alpha=.3,
                     visible = False,linestyle = '-',linewidth = 3,
                     label = r'$\tau$'+ ' = Turn \n'+r'$\Omega$'+' = Other'))
    ax.legend(legs,
              map(lambda x:x[0].get_label(),legs)
              ,labelspacing =2,
              numpoints = 1)        

            #ax.arrow(x,y,dx,dy, head_width = hwidth, width = hwidth/5,alpha = aalpha,zorder = -5)
def nussinov(sequence):
    #pairing matrix
    'G,C,A,U'
    M=array([[0,1,0,1],
             [1,0,0,0],
             [0,0,0,1],
             [1,0,1,0]])
    skey = {'G':0,
            'C':1,
            'A':2,
            'U':3}
    l =len(sequence)
    S = zeros([l,l])
    TB = zeros([l,l],int) -4

    for i in range(l)[::-1]:
        for j in range(i+2,l):
            e0 = skey[sequence[i]]
            e1 = skey[sequence[j]]

            max_k = -1
            max_kscore = 0
            for k in range(i+1,j-1):
                score = S[i,k]+S[k+1,j]
                if score > max_kscore or max_k ==-1:
                    max_kscore = score
                    max_k = k
                    
            cell_scores=[ S[i+1,j-1]+M[e0,e1],
                          S[i+1,j],
                          S[i,j-1],
                          max_kscore]
            best_cell = argmax(cell_scores)
            score = cell_scores[best_cell]

            S[i,j] = score
            TB[i,j] = (-1*best_cell)
            if best_cell == 3:
                TB[i,j] = max_k
            
    return S, TB

def nussinov_traceback(sequence, M, TB):
    l = len(sequence)
    i = 0
    j = l - 1

    ptrs =[(i,j,[0])]
    pairs = []
    bif_count = 1
    while ptrs:
        i,j,bif = ptrs.pop()
        t = TB[i,j]
        if t ==0:
            ptrs.insert(0,(i+1,j-1,bif))
            pairs.append((i,j,bif))
        elif t==-1:
            ptrs.insert(0,(i+1,j,bif))
        elif t==-2:
            ptrs.insert(0,(i,j-1,bif))
        elif t==-4:
            #we're in the middle, continue.
            continue
        else:
            #bifuricate
            
            bif0 = list(bif)
            bif0.append(bif_count+1)
            bif1 = list(bif)
            bif1.append(bif_count+2)

            ptrs.append((i,t,bif0))
            ptrs.append((t+1,j,bif1))
            bif_count += 2
        

    return pairs
def draw_rna(seq,pairs_unsorted):
    f = plt.figure(2)
    f.clear()
    rect = [-1,1,-1,1]
    ax = f.add_axes([0,0,1,1],xlim=rect[0:2],ylim=rect[2:4],
                    aspect='equal')
    l = len(seq)
    xs,ys = [],[]

    pairs =sorted(pairs_unsorted,
                  key = lambda x:x[2][-1])
    max_bif = max(map(lambda x:x[2][-1],pairs))
    bif_colors = getct(max_bif+1)
    bif_colors = [[0,0,1]] *(max_bif+1)
    from matplotlib.patches import Circle
    cir = Circle((0,0),radius=.5,
                 facecolor = '1',edgecolor='.8',
                 zorder = -1)
    ax.add_patch(cir)


    pairs =sorted(pairs_unsorted,
                  key = lambda x:x[2][-1])
    for i in arange(float(l)):
        theta = i*2*pi/l
        r=.5
        x = r*cos(theta)
        y= r*sin(theta)
        xs.append(x)
        ys.append(y)

        #ax.scatter(x,y,60,
        #           facecolor ='0',
        #           edgecolor ='.5',
        #           alpha = 1)
    for p in pairs:
        t0 = mod_with_min((p[0]) * 2*pi / l,2*pi,-pi)
        t1 = mod_with_min((p[1]) * 2*pi / l,2*pi,t0)

        rtval = mod_with_min(t1-t0,  2*pi, -pi)
        rtval0 = sign(rtval)*(pi - abs(rtval))/pi

        cstr = "arc3,rad="+str(rtval0)
        ax.annotate('',
                    xy=[xs[p[0]],ys[p[0]]],xycoords='data',
                    xytext=[xs[p[1]],ys[p[1]]],textcoords='data',
                    
                    arrowprops=dict(arrowstyle="wedge",
                                    shrinkA = 10,
                                    shrinkB = 5,
                                    fc=bif_colors[p[2][-1]], ec="none",
                                    connectionstyle=cstr)
                    )
        rtval1 = sign(rtval)/pi *pow(abs(rtval) + pi,1)
        
    ax2 = f.add_axes([0,0,1,1],frame_on=False,polar = True)
    ax2.set_rlim([0,1.0])
    ax2.set_axis_off()
    grpcounts = {}

    extents = zeros(max_bif+1)
    hashes = zeros([max_bif+1,l],float)
    xax = arange(float(l))/l*2 * pi
    for k,grp in it.groupby(pairs, lambda x:x[2][-1]):
        glist = list(grp)
        n = len(glist)
        grpcounts[k] = n
        sofar = 0
        for p in glist:
            
            p01 = array(sorted(p[0:2]))
            print 'p0'+str(p01)
            if mean(abs(mod_with_min(p01,l,-l/2))) < l/4 and p01[1]-p01[0] > 50:
                incl0 = True
            else:
                incl0 = False

            if incl0:
                elts = mod_with_min(arange(p01[1],p01[0]+l),l,0)
            else:
                elts = arange(p01[0],p01[1])
    
            n = len(elts)
            
            trippy = False
            decay = 1
            do_decay = True
            if not do_decay:
                xvals = zeros(n)+1
            elif trippy:
                xvals = 1/(1+exp(linspace(-1,1,n)))
            else:
                xvals = exp(-(1-abs(linspace(-1,1,n)))*decay)
            hashes[int(k),elts] += xvals
            extent = len(elts)
            if extent > extents[k]:
                extents[k] = extent

            t0 = mod_with_min((p[0]) * 2*pi / l,2*pi,-pi)
            t1 = mod_with_min((p[1]) * 2*pi / l,2*pi,t0)

            arc_angles = array([t0,t1])*180/pi
            if t1 - t0 > pi: arc_angles = arc_angles[::-1]
        

            from matplotlib.patches import Arc
            arc = Arc([0,0],1.1,1.1,0,
                      arc_angles[0],arc_angles[1],
                      zorder = len(p[2]),
                      color = bif_colors[p[2][-1]],
                      linewidth =len(p[2])*7,
                      )
            #ax2.add_patch(arc)
            sofar+=1
    hsum = zeros(l)
    
    hashkeys = argsort(extents)[::-1]
    max_sum = max(np.sum(hashes,0))
    inds = mod(arange(l+1),l)
    for k in hashkeys:
        h = hashes[k]
        spacer = greater(h,0)*.05
        ax2.fill_between(xax[inds],((hsum+h)/max_sum/2 +.5)[inds],
                         (hsum/max_sum/2 +.5)[inds],
                         edgecolor = 'none',
                         facecolor = bif_colors[k],
                         alpha = 1)
        
        ax2.fill_between(xax[inds],((hsum+h)/max_sum/2 +.5)[inds],
                         ((hsum+h)/max_sum/2 + spacer+.5)[inds],
                         edgecolor = 'none',
                         facecolor = 'white',
                         alpha = 1)
                
        hsum+=h+spacer
def mod_with_min(x_in,modby,minval):
    x = array(x_in)
    x -= minval
    x = mod(x,modby)
    x += minval
    return x

def draw_sum(M,N,T,E,skey,ekey):
    nt, l = M.shape
    syms = sym_map()

    items = list(skey.itervalues())
    smap = [  syms[skey.keys()[items.index(i)]] for i in range(max(items)+1)  ]
    items = list(ekey.itervalues())
    emap = [  syms[ekey.keys()[items.index(i)]] for i in range(max(items)+1)  ]
    
    print 'drawing'
    f = plt.figure(3)
    f.clear()
    rect = array([-1.0,l,-.5,1.5])
    ax2 = f.add_axes([0,0,1,1],xlim = rect[0:2],ylim=rect[2:4])
    ax2.autoscale(False)
    
    ct = getct(len(smap))
    for i in range(l):
        psum = np.sum(M[:,i])
        for j in range(nt):
            x = i
            y = float(j) / (nt -1)
            sym = smap[j]
            ax2.scatter(x,y,800,color = 'white',
                        edgecolor = ct[j],linewidth=5,
                        alpha = M[j,i]/psum)
            
            ax2.annotate(sym,(x,y),
                         size = 20,
                         ha = 'center',
                         va = 'center')
            ax2.annotate(str(int(N[j,i])) ,(x,y),
                         xytext = (20,20),textcoords = 'offset points')
            
            if i == l -1:
                path_prob=M[j,i]
                ax2.annotate('P = '+str(path_prob), xy = (x,y),xytext =(20,5),  textcoords = 'offset points')
                         
    
    ax2.annotate('''Opacities denote likelihood of any path transitioning through a state at any given time
Numbers denote number of paths transitioning through a state at a given time.''',
                 xy=(.1,.15),size = 16, xycoords = 'axes fraction',va = 'top')
    ax2.annotate('''Total probability: ''' + str(np.sum(M[:,l-1])) + '''
Total paths: ''' + str(sum(N[:,l-1])) + '''
''',
                 xy=(.77,.15),size = 16, xycoords = 'axes fraction',va = 'top')
    ax2.annotate('All Paths', 
                 xy=(.1,.9),size = 25, xycoords = 'axes fraction')
         
    plt.show()
    
def sym_map():
    syms={'M':'M',
          'A':'A',
          'L':'L',
          'N':'N',
          'G':'G',
          'Alpha Helix':r'$\alpha$',
          'Beta Sheet':r'$\beta$',
          'Other':r'$\Omega$',
          'Turn':r'$\tau$',
          'E':'E'
        }
    return syms 
def draw_path(path,sequence,M,T,E,skey,ekey):
    l = len(sequence)
    nt = shape(T)[0]
    syms = sym_map()

    items = list(skey.itervalues())
    smap = [  syms[skey.keys()[items.index(i)]] for i in range(max(items)+1)  ]
    items = list(ekey.itervalues())
    emap = [  syms[ekey.keys()[items.index(i)]] for i in range(max(items)+1)  ]

    f = plt.figure(0)
    f.clear()
    rect = array([-1.0,l,-.5,1.5])
    ax2 = f.add_axes([0,0,1,1],xlim = rect[0:2],ylim=rect[2:4])
    ax2.autoscale(False)
    
    ax2.annotate('''Opacities denote likelihood of the best path to time 't' transitioning through a state at any given time.
(i.e, the value of the Viterby table)''',
                 xy=(.1,.1),size = 16, xycoords = 'axes fraction')
    ax2.annotate('Viterby path', 
                 xy=(.1,.9),size = 25, xycoords = 'axes fraction')



    ct = getct(len(smap))
    for i in range(l):
        psum = np.sum(M[:,i])
        for j in range(nt):
            sym = smap[j]
            y_chosen =  smap[int(path[i])] == sym 
            print y_chosen

            x = i
            y = float(j) / (nt -1)

            if y_chosen:
                ax2.scatter(x,y,1200,color = 'none',
                            edgecolor = '0',linewidth=1,
                            alpha = 1)

                path_prob=M[j,i]
                ax2.annotate('P = '+str(path_prob), xy = (x,y),xytext =(30,10),  textcoords = 'offset points')
            ax2.scatter(x,y,800,color = 'white',
                        edgecolor = ct[j],linewidth=5,
                        alpha = M[j,i]/psum)



            ax2.annotate(sym,(x,y),
                         size = 20,
                         ha = 'center',
                         va = 'center')
            if i >0 and y_chosen:
                ax2.annotate('',(x,y),xycoords = 'data',
                             xytext = (xl,yl),textcoords='data',
                             arrowprops=dict(arrowstyle='wedge',
                                             shrinkA = 50,
                                             shrinkB = 50,
                                             fc = 'none'),
                             zorder = -1
                             )
            if y_chosen :
                xl = x
                yl = y
def init_p2():
    M1 = array([[.180,.274,.426,.120],
          [.171,.367,.274,.188],
          [.161,.339,.375,.125],
          [.079,.355,.384,.182]])
    M2 = array([[.3,.205,.285,.210],
          [.322,.298,.078,.302],
          [.248,.246,.298,.208],
          [.177,.239,.292,.292]])

    
    return M1,M2
def p2_randomseq(M, l = 6, seed = 2):
    n = 4
    Mcum = zeros((n,n))
    for i in range(n):
        for j in range(n):
            Mcum[i,j] = np.sum(M[i,:j])
            
    random.seed(seed)
    u = random.uniform(0,1)
    cur = int(floor(u/.25))
    out = [cur]
    for i  in range(l -1):
        cum = Mcum[cur]
        u = random.uniform(0,1)
        
        cur = nonzero(less_equal(cum,u))[0][-1]
        out.append(cur)
    out = array(out)
    return out
def p2_prob(M,seq):
    l = len(seq)
    prob = np.product(M[seq[:-2],seq[1:-1]])
    return prob
    
def p2_drawprobs( gen1, gen2 ):
    
    maxprob = np.max([gen1[0],gen1[1],gen2[0],gen2[1]])
    gen1/= maxprob
    gen2/= maxprob

    f = plt.figure(0)
    f.clear()
    ax1 = f.add_axes([.05,.05,.9,.4],xlim = [0,1])
    ax2 = f.add_axes([.05,.55,.9,.4],xlim = [0,1])

    n = len(gen1[0])
    xax = linspace(.05,.95,n)

    for i in range(2):
        if i ==1:
            ax = ax1
            correct = gen1[0]
            wrong = gen1[1]
            name = 'CpG sequences'

        else:
            ax = ax2
            correct = gen2[1]
            wrong = gen2[0]
            name = 'Non CpG sequences'


        args = argsort(  correct)[::-1]
        ax.fill_between(xax,correct[args],wrong[args],
                        interpolate = True,
                        where = greater(correct[args],wrong[args]),color = 'blue')
        ax.fill_between(xax,wrong[args],correct[args],
                        interpolate = True,
                        where = greater(wrong[args],correct[args]),color = 'red')
        ax.plot(xax,wrong[args],linewidth =3,color = '1')
        cplot = ax.plot(xax,wrong[args],linewidth = 1,linestyle = '--',color = '0')

        ax.plot(xax,correct[args],linewidth =11,color = '1')
        wplot = ax.plot(xax,correct[args],linewidth = 5,color = '0')
        
        ax.annotate(name,xy = [.1,.8],xycoords = 'axes fraction', size = 24)
        ax.annotate('''Classification error: '''+str(float(len(nonzero(greater(wrong,correct))[0]))/len(wrong))+'''
Mean classification certainty: ''' + str(mean(( correct - wrong)[nonzero(greater(correct,wrong))[0]])), xy = (.5,.8),
                    ha = 'center',
                    va = 'top',
                    xycoords = 'axes fraction')
        ax.legend((cplot,wplot),
                 ('Correct classifier probability',
                  'Incorrect classifier probability'))
        
def p4_plot(scores):
    xs = scores['l']
    ys = scores['scores']


    print xs
    print ys
    f = plt.figure(0)
    f.clear()
    ax = f.add_axes([0,0,1,1])
    ax.plot(xs,ys)
