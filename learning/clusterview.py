import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import utils.colors as mycolors
import netutils as nu


def specview(means, membership, fig = 12):

    f = plt.figure(fig)
    f.clear()
    ax =  f.add_axes([.05,.05,.9,.9])
    
    #ax.plot(membership[np.argsort(membership)])

    img = means
    activity = mycolors.imgcolor(img,BW = True)
    ax.imshow(img,aspect = 'auto')

def viewmany(all_means, all_clusters, fig = 12):
    n = len(all_means)
    f = plt.figure(fig)
    f.clear()
    print '''Running viewmany.py

For now, viewmany assumes that k is equal across clustering instances
this is not really important but has to do with how TF projections are
stored.
'''
    #1 k.
    k = len(all_means[0])

    ax1 = f.add_axes([.05,.05,.95,.4])
    ax2 = f.add_axes([.05,.55,.95,.4])
    ct0 = mycolors.getct(n)

    sqa = nu.net_square_affinity()[0]
    aff = nu.net_affinity()[0]

    #tf_sqidxs should have length = ntf
    #with each element giving the coordinate of the
    #i'th tf in sqa space.


    sqidxs = nu.net_sq_keyidxs()
    n_tfidxs = nu.net_tf_keyidxs() 
    trgs,tfs = nu.parse_net()
    tf_sqidxs = [sqidxs[key] for key in tfs.keys()]
    tfidxs = n_tfidxs.values()
    ntf = len(tfidxs)

    tfweights = zeros(ntf,int)
    #find tfs of general interest, choosing at most ten for each clustering
    ntf_each = 20
    
    print '''...Computing representative TFs for each clustering.

In the current formulation, we project each mean on to associated tf
and then normalize each projection so that each mean has equal weight
in TF selection.

Not that we have handled the case where we have clusted in TF space
explicitly (e.g, dim = 541) and where we are in gene space explicitly,
(e.g., dim = 8321, GG matrix or svdU). svdV is emphatically not handled.
Neither would svdU of TF-TF which is actually the the exact same thing.'''
    

    
    TFprojs= zeros((n,k,ntf))
    for i in range(n):
        m = all_means[i]
        dim = shape(m)[1]
        #we are now going to project clusters on to the tfs
        #in this form, we only need rows corresponding to tfs.

        if dim> 500:
            #If dim = 541, we just read off the most important tfs
            this_tf_sum = np.abs(m[:,tfidxs])
            TFprojs[i,:,:] = this_tf_sum
            #normalize clusters
            this_tf_sum = this_tf_sum / np.sum(this_tf_sum,1)[:,newaxis]
            this_tf_sum = np.sum(this_tf_sum,0)
    
        #Now, since we are at the moment only working with GG
        #and SVD_U, we are in gene space and can undo the mapping
        #with sqaT
        elif dim > 8000:
            #remember, ROWS of the matrix correspond to the
            #target space.
            a = sqa.T[tf_sqidxs,:]            
            this_tf_sum = np.abs(np.sum(a[newaxis,:,:]*m[:,newaxis,:],2))
            TFprojs[i,:,:] = this_tf_sum
            #normalize so that each mean has the same weight
            this_tf_sum = this_tf_sum / np.sum(this_tf_sum,1)[:,newaxis]
            #sum over cluster means to find the most important tfs
            this_tf_sum = np.sum(this_tf_sum,0)
            
    

        best = argsort(this_tf_sum)[::-1]
        tfweights[best[0:ntf_each]]=1
    print '''Finished computing representative TFs
'''

    tfs_of_interest = nonzero(tfweights)[0]
    ntf = len(tfs_of_interest)
    avg_unshared = float(ntf)/(n * ntf_each)
    avg_shared = 1. - float(ntf)/(n * ntf_each)
    print '''Allowing for each cluster to choose '+str(ntf_each) + 'tfs,
we got ''' + str(ntf) + ''' tfs of interest.
or a mean sharing ratio of ''' + str(round(avg_shared,3))+ '''.'''

    #get a color table for clusters.
    ct = mycolors.getct(n)

    for i in range(n):
        #p stands for 'point' as in datapoint.
        #data points are labeled with clusters.

        xax = linspace(0,1,ntf)

        ax1.plot(xax,np.sum(TFprojs[i,:,tfs_of_interest],1)/np.max(TFprojs[i,:,tfs_of_interest],1),color = ct[i])

    return TFprojs




def one(all_means, all_mems,
        tfp, 
        axis = 'tf',
        idxs = [0,1], fig = 5
        ,choice_ax = 'x'
        ,nrml = 'axis'
        ,sorting = 'axis'):
    m = all_means[idxs[0]]
    c = all_mems[idxs[0]]
    proj=abs(tfp[idxs[0],:,:])

    m2 = all_means[idxs[1]]
    c2 = all_mems[idxs[1]]
    proj2=abs(tfp[idxs[1],:,:])

    
    sqidxs = nu.net_sq_keyidxs()
    n_tfidxs = nu.net_tf_keyidxs() 
    trgs,tfs = nu.parse_net()
    tf_sqidxs = [sqidxs[key] for key in tfs.keys()]
    gene_sqidxs = [sqidxs[key] for key in trgs.keys()]

    tfk = nu.net_tf_keyidxs()
    tgk = nu.net_trg_keyidxs()
    tf_aidx = [ tfk[key] for key in tfs.keys()]
    gene_aidx = [ tgk[key] for key in trgs.keys()]


    tfidxs = tf_aidx
 
    k = len(m)
    ntf = len(tf_sqidxs)
    ng = len(gene_sqidxs)


    print '''Getting ready to plot clusters mapped on to tf components.

--note--
In its current incarnation, netutils orders tfs by their out degree
and genes by their in degree.

Thus viewmany() orders projects by TF out degree. Left unsorted, this
is the order of the TF x axis.'''
    
    #how to normalize the image?
    #axis: equal sum for each tf over all clusters.
    #other: equal sums for each cluster in img


    
    nrml = 'axis'
    nrml_type = lambda x,y:np.max(x,y)

    sorting = 'other'
    
    
    print axis

    d0 = shape(m)[1]
    d2 = shape(m2)[1]

    show_membership = True

    if axis == 'tf':    

        if sorting == 'axis':
            img = proj
            mean_tfval = argmax(img,1) 
            c_srt = np.argsort( mean_tfval)


            img = img[c_srt,:]    

            img2 = proj2
            mean_tfval = argmax(img2,1) 
            c_srt = np.argsort( mean_tfval)

            img2 = img2[c_srt,:]    
        else:
            img = proj
            mean_tfval = argmax(img,0) 
            c_srt = np.argsort( mean_tfval)


            img = img[:,c_srt]    

            img2 = proj2
            mean_tfval = argmax(img2,0) 
            c_srt = np.argsort( mean_tfval)

            img2 = img2[:,c_srt]    
    elif axis =='gene':  
        maxgene = 200
        gsort = argsort(c)
        
        
        if d0 == 8321 and not show_membership:
            img = m[:,gsort][:,:maxgene]
        else:
            img = zeros((k,ng))
            for i in range(ng):
                img[c[i],i] = 1
        if d2 == 8321 and not show_membership:

            img2 = m2[:,gsort][:,:maxgene] 
        else:
            img2 = zeros((k,ng))
            for i in range(ng):
                img2[c2[i],i] = 1


    #normalize to generate an image
    if nrml == 'axis':
        img2 = img2/nrml_type(img2,0)[newaxis,:]
        img = img/nrml_type(img,0)[newaxis,:]
    else:
        img2 = img2/nrml_type(img2,1)[:,newaxis]
        img = img/nrml_type(img,1)[:,newaxis]       
            
    
    img /= np.max(img)
    img2 /=np.max(img)

    img_show= img[:,:,newaxis] *[0,0,1] + img2[:,:,newaxis]*[1,0,0]
    


    f = plt.figure(fig)
    f.clear()


    ax = f.add_axes([.05,.05,.9,.9])
    ax.imshow(img_show[:,:,:], aspect = 'auto')

    nc = shape(img)[0]
    xs, ys, rs, cs = [[] for i in range(4)]
    
    nchoice = 1
    if choice_ax == 'y':

        dim =  shape(img)[0]
        maxes = [argsort(img,1)[::-1][:,:nchoice],
                 argsort(img2,1)[::-1][:,:nchoice]]
    elif choice_ax == 'x':
                     
        dim =  shape(img)[1]
        maxes = [argsort(img,0)[::-1][:nchoice,:],
                 argsort(img2,0)[::-1][:nchoice,:] ]    
    else:
        raise Exception('bad axis')                

    ct = mycolors.getct(len(maxes))
    for j in range(len(maxes)):
        for i in range(dim):
            for k in range(nchoice):
                if choice_ax == 'x':
                    ys.append(maxes[j][k][i])
                    xs.append(i)
                elif choice_ax =='y':
                    xs.append(maxes[j][i][k])
                    ys.append(i)
                else:
                    raise Exception('bad axis')

                rs.append(20 + 30*(1-j))
                cs.append(ct[j])
            
    xs, ys, rs, cs  = np.array(xs),np.array(ys),np.array(rs),np.array(cs)

    ax.scatter(xs,ys,200,'1',edgecolor = 'none')
    ax.scatter(xs,ys,rs,cs,alpha = .8, edgecolor = 'none')
        

    
