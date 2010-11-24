def run(part = '4a'):
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
            
                
def p4_plot(scores):
    xs = scores['l']
    ys = scores['scores']


    print xs
    print ys
    f = plt.figure(0)
    f.clear()
    ax = f.add_axes([0,0,1,1])
    ax.plot(xs,ys)

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
