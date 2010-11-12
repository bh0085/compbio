import pickle
import numpy as np
import matplotlib.pyplot as plt

#Assume that relevant data has been saved in
#the working directory for each problem in the
#form probX.pickle

def prob2():
    import ps1_dotplot_modified as dp
#    dp.main(plotfile = '2a.jpg',
#            file1 = 'data/human-hoxa-region.fa', 
#            file2 = 'data/mouse-hoxa-region.fa',
#            skip = 1,
#            delskip = 0,
#            kmerlen = 30,
#            allow_inversions = False)
#    dp.main(plotfile = '2bi.jpg',
#            file1 = 'data/human-hoxa-region.fa', 
#            file2 = 'data/mouse-hoxa-region.fa',
#            skip = 1,
#            delskip = 0,
#            kmerlen = 100,
#            allow_inversions = False)
#    dp.main(plotfile = '2bii.jpg',
#            file1 = 'data/human-hoxa-region.fa', 
#            file2 = 'data/mouse-hoxa-region.fa',
#            skip = 2,
#            delskip = 0,
#            kmerlen = 60,
#            allow_inversions = False)
#    dp.main(plotfile = '2biii.jpg',
#            file1 = 'data/human-hoxa-region.fa', 
#            file2 = 'data/mouse-hoxa-region.fa',
#            skip = 3,
#            delskip = 0,
#            kmerlen = 90,
#            allow_inversions = False)
#    dp.main(plotfile = '2biv.jpg',
#            file1 = 'data/human-hoxa-region.fa', 
#            file2 = 'data/mouse-hoxa-region.fa',
#            skip = 4,
#            delskip = 0,
#            kmerlen = 120,
#            allow_inversions = False)
#    dp.main(plotfile = '2bv.jpg',
#            file1 = 'data/human-hoxa-region.fa', 
#            file2 = 'data/mouse-hoxa-region.fa',
#            skip = 1,
#            delskip = 3,
#            kmerlen = 100,
#            allow_inversions = False)
#

    dp.main(plotfile = '2e.jpg',
            file1 = 'data/human-hoxa-region.fa', 
            file2 = 'data/human-hoxa-region-modified.fa',
            skip = 1,
            delskip = 0,
            kmerlen = 30,
            allow_inversions = True)
    


def prob3():
    data = pickle.load(open('prob3.pickle'))
    hexamers = data['hexamer']
    nums = np.array(data['nums'])
    cons = np.array(data['cons'])
    frac = np.array(data['frac'])
    pvals = np.array(data['pvals'])

    fig1 = plt.figure(1,[7,4])
    fig1.clf()
    ax0 = plt.subplot(2,2,1)
    ax1 = plt.subplot(2,2,2)
    ax2 = plt.subplot(2,2,3)
    ax3 = plt.subplot(2,2,4)
    ax0.set_title('pvals')
    ax1.set_title('fractional conservation sorted by pval')


    ax0.set_yscale('log')
    ax0.set_ylim(1e-8,1)
    ax1.set_yscale('log')

    from mymath import ssmooth as ss

    sp = np.argsort(pvals)
    sf = np.argsort(frac)
    sn = np.argsort(nums)
    ax0.plot(pvals[sp])
    #ax1.plot(ss.smooth(frac[sp],10))
    ax1.plot(frac[sp])

    ax2.plot(pvals[sp], ss.smooth(nums[sp],5))
    ax3.plot(frac[sf], ss.smooth(nums[sf],40))

    n = len(pvals)
    sig = .05
    sig_b = .05 / 4096
    where_sig = np.nonzero(pvals < sig)[0]
    where_sig_b = np.nonzero(pvals < sig_b)[0]
    
    nsig = len(where_sig)
    nbon = len(where_sig_b)

    from matplotlib.backends.backend_pdf import PdfPages
    pp  = PdfPages('p3a.pdf')
    pp.savefig(ax0.figure)
    pp.close()

    #Make a pie chart of the number of elements at a given significance

    plt.figure(2,figsize = [5,5])
    plt.clf()
    piechart = plt.pie([nbon, nsig - nbon, n - nsig - nbon],
                       explode = np.zeros(3)+.03,
                       labels = ['B. sig:\n '+str(nbon),
                                 'sig (p = .05):\n '+str(nsig),
                                 'Other: '+str(n - nsig - nbon)])

    #Draw a histogram plotting relative abundances of AT and C/G in
    #significant and/or highly conserved motifs.
    atfrac = []
    at = 0.0
    for item in where_sig_b:
        at += hexamers[item].count('A') + hexamers[item].count('T')
    atfrac.append(at / (len(where_sig_b) * 6))
    
    at = 0.0
    for item in where_sig:
        at += hexamers[item].count('A') + hexamers[item].count('T')
    atfrac.append(at / (len(where_sig) * 6))
    
    at = 0.0
    for item in sf[range(-1,-101,-1)]:
        at += hexamers[item].count('A') + hexamers[item].count('T') 
    atfrac.append(at / (100 * 6))
    
    at = 0.0
    for item in sn[range(-1,-101,-1)]:
        at += hexamers[item].count('A') + hexamers[item].count('T')
    atfrac.append(at / (100 * 6))
    at = 0.0
    at_all = []
    colors = []
    for item in range(len(hexamers)):
        atcount =  hexamers[item].count('A') + hexamers[item].count('T')
        at += atcount
        at_all.append(atcount)
        colors.append([float(atcount)/6, 0,0])
    mean_at = at/len(hexamers)/6

    plt.figure(3,  [5,3])
    plt.clf()
    ax = plt.subplot(1,1,1)
    ax.plot(range(4),atfrac, 'r')
    ax.plot(range(4),np.zeros(4) + mean_at, 'b')

    plt.figure(4,[5,5])
    plt.clf()
    ax = plt.subplot(1,1,1)
    ax.set_ylim([.5*1e-10,1])
    ax.set_xlim([5e1,5e4])
    ax.set_xscale('log')
    ax.set_yscale('log')


    nidxs = sf[range(-1,-101,-1)]
    pvals_rounded = np.array(map(lambda x: (x> 1 and 1) or x > 1e-10 and x or x <  1e-10 and 1e-10,pvals))

    ax.fill_between([1,1e6],np.zeros(2) + .05,np.zeros(2) + 1e-11 ,alpha = 1,
                    edgecolor = [.3,.3,.3],facecolor =  [.95,.95,.95],linestyle = 'solid',linewidth = 1)

    ax.fill_between([1,1e6],np.zeros(2) + .05/4096,np.zeros(2) + 1e-11 ,alpha = 1,
                    edgecolor = [.2,.2,.2],facecolor =  [.8,.8,.8], linestyle = 'solid',linewidth = 1)


    ax.scatter(nums, pvals_rounded, s =400, c = 'w',alpha = .2,edgecolor = [0,0,0],facecolor = 'none')    
    ax.scatter(nums, pvals_rounded, s =400, c = 'w',edgecolor = 'none')

    ax.scatter(nums[nidxs], pvals_rounded[nidxs], s= 100,c='b',alpha = .2)

    ax.scatter(nums, pvals_rounded, s =30, c = colors, edgecolor = 'none')
    ax.set_xlabel('Number of occurences of Hexamer')
    ax.set_ylabel('Conservation p-value')

    legs = (plt.plot(1,1,markeredgecolor = 'none',markerfacecolor = [1,0,0],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'hexhigh'),
            plt.plot(1,1,markeredgecolor ='none',markerfacecolor = [0,0,0],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'hexhigh'),
            plt.plot(1,1,markeredgecolor = 'none',markerfacecolor = [.6,.6,1],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'hexhigh'))

    ax.legend(legs,
              ('High AT Hexamer',
               'Low AT Hexamer',
               'Top 100 Conserved')
              ,numpoints = 1,markerscale = 2)
    ax.annotate('p < .05',[1e4,1e-3] )
    ax.annotate('5% Confidence with Bonferroni\n(p < .05/4096)',[5e3,1e-7])
    
    wat = np.nonzero(map(lambda x : x == 'ATATAT' and 1 , hexamers))[0][0]
    wtt = np.nonzero(map(lambda x : x == 'TTTTTT' and 1 , hexamers))[0][0]
    wcc = np.nonzero(map(lambda x : x == 'CCGGGT' and 1 , hexamers))[0][0]
    wcg = np.nonzero(map(lambda x : x == 'CGCGCG' and 1 , hexamers))[0][0]
    ax.annotate('ATATAT',[nums[wat],pvals_rounded[wat]],
                color= [0,0,0], bbox=dict(boxstyle="round", fc=[1,1,1]))
    ax.annotate('TTTTTT',[nums[wtt],pvals_rounded[wtt]],
                color= [0,0,0], bbox=dict(boxstyle="round", fc=[1,1,1]))

    ax.annotate('CCGGGT',[nums[wcc],pvals_rounded[wcc]],
                color= [0,0,0], bbox=dict(boxstyle="round", fc=[1,1,1]))

    ax.annotate('CGCGCG',[nums[wcg],pvals_rounded[wcg]],
                color= [0,0,0], bbox=dict(boxstyle="round", fc=[1,1,1]))
    
    #raise Exception
    
    #raise Exception
