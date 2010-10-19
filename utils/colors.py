import random

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

def pmcolor(val,threshold = 0):
    if val < -1*threshold:
        return [1,0,0]
    elif val > threshold:
        return [0,0,1]
    else:
        return [0,0,0]
