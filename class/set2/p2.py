from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import random
def run(part = '2a'):
    if part =='2a':
        M1,M2 = init_p2()
        seqs1, seqs2 = [], []
	#set n as you please, right now n = 200.
        for i in range(200):
            seqs1.append( p2_randomseq(M1,6,i))
            seqs2.append( p2_randomseq(M2,6,i+10000))
            
        
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
        
