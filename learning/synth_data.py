from numpy import *
def synth_data(ntg, max_unshared, ntf_shared):
        
    unshared_per_tg = random.random_integers(0,max_unshared-1,ntg)
    all_shared_tfs = range(ntf_shared)
    
    #30 training points
    nt = 30
    #synthetic tf expr values
    x_expr= list(random.uniform(-1,1,(ntf_shared,nt)) )
    null = ( [x_expr.extend(list((random.uniform(-1,1,(max_unshared, nt)))))    \
                        for i in range(ntg)    ])
    x_expr = array(x_expr)
    
    cxns = []
    unshared_weights = []
    for i in range(ntg):
        #connect to a random set of unshared tfs
        nu = unshared_per_tg[i]
        #connect to all shared tfs
        ns = ntf_shared

        which_shared = range(ns)
        which_unshared = random.random_integers(0,max_unshared-1 , nu)
        these_s_cxns = which_shared  
        these_u_cxns = [ x + ntf_shared + max_unshared * i for x in which_unshared]
        #create a distinct unshared weight for each tg/tf
        unshared_weights.append(random.uniform(-1,1,unshared_per_tg[i])) 

        cxns.append([these_s_cxns,these_u_cxns])

    #create shared weights for shared tfs
    shared_weights = random.uniform(-1,1,ntf_shared)

    y_expr = []
    for y in range(ntg):
        these_weights =array( list(shared_weights) + list(unshared_weights[y]))
        these_cxns = cxns[y][0] + cxns[y][1]
        these_wts = array(list(shared_weights[:]) +list( unshared_weights[y]))

        this_expr = sum(x_expr[these_cxns,:]*these_wts[:,newaxis],0)
        y_expr.append(this_expr)
        
    #synthetic y_expr values, [tf] * [weights].
    y_expr = array(y_expr)

    return x_expr, y_expr


