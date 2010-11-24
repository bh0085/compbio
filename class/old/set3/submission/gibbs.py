#!/usr/bin/env python
import sys
import string
import random

import matplotlib.pyplot as plt
from numpy import *
import numpy as np

#### INSTRUCTIONS FOR USE:
# Don't call from the command line... call run or runall from within python.


alphabet = ['A', 'G', 'C', 'T']
abcdict = {'A':0,
           'G':1,
           'C':2,
           'T':3}

#Runs gibbs sampling to a max number of iterations given by max_r
#Many free parameters... Check inside.
#
#Hardcoded parameters choose things like whether to plot a score.
#
def GibbsSampler(S0,L,  max_r = 1000):
    ns = len(S0)
    na = len(alphabet)
    random.seed()
    
    #convert each sequence to an [n x 4] array giving a positional representation
    #of nucleotide content
    S = []
    for seq in S0:
        n = len(seq)
        arr = zeros((n, na),int)
        for i in range(n):
            arr[i][abcdict[seq[i]]] = 1
        S.append(arr)

    #Get a random list of subsequences to start.
    subs = zeros((ns, L, na))
    for i in range(ns):
        start = random.random_integers(0,len(S[i]) - L -1)
        subs[i] = S[i][range(start,start+L), :]
      

    
    #PARAMETERS:
    #Plot?
    do_plot = False
        
    #randomize subsequence omit choice.
    omit_random = False
    #or increment an omission pointer.
    last_omit = 0

    #Simulated annealing adds a bunch of noise that decreases
    #as the simulation progresses.
    #It seems to work!
    sim_anneal = True
    t_init = .2
    #Just noise, no simulated annealing.
    pwm_noise = False

    #weight motifs is slower to get new ideas.
    #its a bad idea.
    weight_motifs = False

    #PWM makes the pwm out of the "mode"
    #I guess this is another thing that doesn't really work.
    pwm_mode = False

    #Score Selection:
    #just_max seems to pretty much break the algorithm. 
    #too easy to get stuck
    just_max = False
    


    converged = 0
    runs = 0 
    max_scores = []  
    
    #Convergence testing
    cvg_window = ns* 2
    cvg_failcount = 0
    cvg_threshold = .05
    cvg_maxfail = 100
    cvg_testmod = 1

    if do_plot:
        f = plt.figure(1)
        f.clear()
        ax1 = f.add_subplot(210)
        ax2 = f.add_subplot(211)

    while not converged:
        if omit_random:
            omit = random.random_integers(0,ns-1)
        else:
            last_omit = mod(last_omit +1, ns)
            omit = last_omit

        pwm = get_pwm(subs,
                      omit=omit,
                      pwm_noise = pwm_noise,
                      sim_anneal = sim_anneal,
                      t_init = t_init,
                      pwm_mode = pwm_mode,
                      runs = runs,
                      max_r = max_r)

        #Get score for subsequence chosen.
        SUB = S[omit]
        scores, cums = score_seq(SUB, pwm)

        #Choose a random sample.
        sampled = choose_random(scores,cums, just_max = just_max)
        subs[omit,:,:] = S[omit][sampled:sampled+L,:] 
        score = scores[sampled]

        #Weight new sequences based on how well they match?
        if weight_motifs:
            subs[omit,:,:] *= (.01 + scores[sampled])
        max_scores.append(scores[sampled])
        
        
        
        if do_plot:
            if mod(runs, max_r / 23)  == 0:
            #this_scr_sort = log(sorted(scores + .0001)[::-1][0:100])
            #ax1.fill_between(range(100),this_scr_sort + last_scr_sort,last_scr_sort,
                            # facecolor = random.uniform(0,1,3))
            #last_scr_sort += this_scr_sort
                ax1.plot((scores))
                ax1.scatter(argmax(scores), (max(scores)),90,color = 'blue')
                ax1.scatter(sampled,score,500*(score + .1),color = [float(runs)/max_r,0,0])
                print argmax(scores), sampled

        runs += 1        

        if mod(runs, cvg_testmod) == 0 :
            if cvg_window*6 < runs: 
                old_win = scores[:-cvg_window*3:-1]
                new_win = scores[:-cvg_window:-1]
                if (mean(new_win) - mean(old_win)) / (mean(old_win)) < cvg_threshold:
                    cvg_failcount += 1
                else:
                    cvg_failcount = 0 

            if cvg_failcount > cvg_maxfail:
                break

        if runs > max_r:
            converged = 1

    pwm = array(sum(subs,0),float)
    pwm/= sum(pwm,1)[:, newaxis]

    if do_plot: ax2.plot(max_scores)
    return pwm, scores





#Choose a random motif given a list of scores and a cumulative distribution
#of scores.
def choose_random(scores, cums, just_max = False):
    if not just_max:
        sampled = nonzero(greater_equal(cums, random.uniform(0,cums[-1])))[0][0]
    else:
        sampled = argmax(scores)
    return sampled

#Get a position weight matrix given a sequence to omit.
#Take noise parameters.
#
#Most interesting noise parameter is sim_anneal which 
#adds a simulation annealing step to weight matrix construction.
#Simulated annealing initial temperature is given by t_init.
def get_pwm(SUBS,
            omit = None,
            pwm_noise = None,
            sim_anneal = None,
            t_init = None,
            pwm_mode = None,
            runs = None,
            max_r = None
            ):
    ns, L, na = shape(SUBS)
    s2 = concatenate((SUBS[:omit],SUBS[omit+1:]))


    #in the case of mode, just choose the
    #best representative at each base position
    if not pwm_mode:
        pwm = array(sum(s2,0),float)
        pwm/= sum(pwm,1)[:, newaxis]


    else:
        maxes =  argmax(sum(s2,0),1)
        max_add = .9
        remains = (1.0 - max_add) / na
        pwm = zeros((L, na)) + remains
        for i in range(L):
            
            pwm[i,maxes[i]] += max_add
    
    if pwm_noise == True:
        pwm +=.01
        pwm/= sum(pwm,1)[:, newaxis] 
    
    if sim_anneal ==True:
        noise = t_init * exp(-1*( (float( runs) / max_r)* 10 ))
        pwm += random.uniform(0,noise, shape(pwm))
        pwm/= sum(pwm,1)[:, newaxis]            
    return pwm


#Calculates scores for different positions in a sequence
#given a pwm.
def score_seq(SUB, PWM):
    n = len(SUB)
    L = shape(PWM)[0]
    n_starts = n - L
    
    scores = zeros(n_starts)
    cums = zeros(n_starts)
        #score subsequences of S_c 
        #and draw from a distribution defined by scores.
    for start in range(n_starts):
        dotprod = SUB[start:start+L,:]*PWM
        score =product(sum(dotprod,1)) 
        scores[start] = score
        if start == 0 :
            cums[start] = score
        else:
            cums[start] = score + cums[start - 1]
    return scores, cums


#Find a consensus motif given a set of scored runs of 
#the motif finder.
def find_motifs(pwms, all_scores, num = 0):
    interval = 50
    pwms_split = []
    scores_split = []

    nd = len(pwms)/interval
    for i in range(nd):
        pwms_split.append(array(pwms[interval*i:interval*(i+1)]))
        scores_split.append(array(all_scores[interval*i:interval*(i+1)]))
        

    f = plt.figure(1)
    f.clear()
    ax1 = f.add_subplot(210)
    ax2 = f.add_subplot(221)
    ax3 = f.add_subplot(222)

    i = num
    pwm = pwms_split[i]
    scores = scores_split[i]
    
    srt = argsort(scores)[::-1]
    srtgood = nonzero(greater_equal(scores[srt], max(scores) * .75))[0]
    srtrealgood = nonzero(greater_equal(scores[srt], max(scores) * .92))[0]

    plots = [ax1.plot(scores[srt]) ,
             ax1.plot(scores[srt[srtgood]],color = 'red'),
             ax1.plot(scores[srt[srtrealgood]],color = 'orange',linewidth = 5)]
    ax1.legend(plots, 
               ['convergence scores over ' + str(interval) + 'runs',
                'runs in the top 25% best scoring',
                'runs in the top 90% best scoring'])

    pbest= squeeze(pwm[srt[0],:,:])

    frac = .1
    pavg = squeeze(sum(pwm[srt[srtrealgood],:,:],0))
        

    ax3.imshow(pavg.T[:,:,newaxis]*[1,0,0], aspect = 'auto', interpolation = 'nearest')
    ax2.imshow(pbest.T[:,:,newaxis]*[1,0,0], aspect = 'auto', interpolation = 'nearest')

    ax2.set_title('Highest scoring motif')
    ax2.set_xlabel('motif position')
    ax2.set_yticklabels(alphabet)
    ax2.set_yticks(range(4))
    ax3.set_title('Average of near optimal scoring motifs')
    ax3.set_xlabel('motif position')
    ax3.set_yticklabels(alphabet)
    ax3.set_yticks(range(4))

    plt.savefig('data'+str(num+1)+'_out.pdf', format = 'pdf')

    print "best score found: "
    print pwm[srt[0]]
    print 
    print "difference between best and second best: "
    print pwm[srt[0]] - pwm[srt[1]]
    print "\n"
    consensus = np.median(pwm,0)
        

#Run
def runall():
    names = ['data1','data2','data3','data4']
    count = 50
    max_r = 1000
    L = 10

    all_pwms = []
    all_scores = []
    for n in names:
        for i in range(count):
            print i
            pwm , scores = run(L =L,datafile = n, max_r = max_r) 
            all_pwms.append(pwm)
            all_scores.append(mean(scores))

    return all_pwms, all_scores
def run(L=10 , datafile = 'data1', max_r = 1000):
    S = readdata(datafile)
    P, Scores = GibbsSampler(S,L,max_r = max_r)
    return P , Scores                         
    

def main():
    L = int(sys.argv[1])
    datafile = sys.argv[2]
    S = readdata(datafile)
	
    P = GibbsSampler(S,L)
	
    print "    ", 
    for i in range(L):
        print "%-5d " % (i+1),
    print ""
	
    for j in range(len(alphabet)):
        print " %s " % alphabet[j], 
        for i in range(L):
            print " %5.3f" % P[j][i],
        print ""
	
def readdata(file):
    data = [];
    for line in open(file,'r'):
        data.append(line[0:-1])
    return data

if __name__ == "__main__": main()

