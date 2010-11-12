import netwriter as nw

def parse_name(name):
    if name == 'unsup':netfile = 'patrick/unsup_patrick.txt'
    if name == 'logistic':netfile = 'patrick/logistic_0.6.txt'
    
    f = open(netfile)
    lines = f.readlines()
    
    tfs = {}
    tgs = {}
    
    for l in lines:
        tf = l.split('\t')[0]
        tg = l.split('\t')[1]
        weight = float(l.split('\t')[2])
        
        tflist = tfs.get(tf, [])
        tflist.append((tg,weight))
        tfs[tf] = tflist

        
        tglist = tgs.get(tg, [])
        tglist.append((tf,weight))
        tgs[tg] = tglist
        
    return tfs, tgs


def get_net(name = 'unsup' , reset = 0):
    hardcopy = True

    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)
    except Exception, e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        net = parse_name(name)
        nw.writenet(name, net, hardcopy = hardcopy)
        return net
