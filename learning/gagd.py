'''This simple implementation of a GA/GD hybrid learner works with a population
of genomes realized from a listed set of posssibilities:

To generate the list of possibilities for a given genome configuration, 
run 

gagd=GAGD()
gagd.initNAMEDSimulation()
gagd.CTemplate()

'''

from pyevolve import G1DList
from pyevolve import GSimpleGA

from pybrain.supervised.trainers import BackpropTrainer
from pybrain import TanhLayer
from pybrain.structure import LinearLayer, SigmoidLayer

import matplotlib.pyplot as plt
from numpy import *
import numpy as np

from compbio.learning import synth_data as sd
from compbio.utils import plots as myplots

def test():
	nhc = 2
        ntg = 2
        ntf_s = 2
        max_tfu = 2
        gagd = GAGD(nhc,ntg,ntf_s, [ max_tfu for i in range(ntg) ] )
        xs, ys = sd.synth_data(ntg,max_tfu,ntf_s)
        g, ga = gagd.sample_genome()
        gagd.init_net()
        gagd.make_cxns_from_genome(g)

        net = gagd.mynn.net
        
        f = plt.figure(0)
        f.clear()
        ax = f.add_subplot(121)
        myplots.draw_pb(ax,net)
        myplots.hideaxes(ax)
        myplots.maketitle(ax,'GANN')
        
        gagd.set_data(xs.T,ys.T)
        gagd.set_trainer()
        gagd.train()
	return



def eval_func(chromosome):
        score = 0.0
   # iterate over the chromosome
        for value in chromosome:
            if value==0:
                score += 1
        return score






class GAGD():
  def makeMotifTemplate(self):
		nhc = self.nhc
		ntg = self.ntg
		ntg_s = self.ntg_s
		ntg_u = self.ntu_u
	
	
	
	
  def initMotifSimulation(self,
				nhc,
				ntg,
				ntg_s,
				ntf_u):
		
		self.nhc = nhc #how many hidden clusters?
		self.ntg = ntg #how many tgs?
		self.ntf_s = ntf_s #how many tfs?
		self.ntf_u = ntf_u #how many unshared tfs per tg.


  def sample_genome(self):

        random.seed()
        hcc_allowed = range(0,self.ntf_s)
        hc_allowed = range(0,self.nhc)
        tg_allowed = range(0,self.ntg)
        hgu_allowed= range(0,max(self.ntf_u))
            #---->>>[range(len(self.ntf_u[i])) for i in range(self.ntg)]
        hgs_allowed = range(0,self.ntf_s)


        n_hc_allowed =  [self.ntf_s]
        n_hgu_allowed = [5]
        n_hgs_allowed = [0]
        n_ogc_allowed = [self.ntg]

        ogc_allowed = range(0,n_ogc_allowed[0])
        ogcc_allowed = range(0,self.ntg)

        w_allowed = [-1.0,0.,1.0]
        ga = []
        g = []

        
        randwt = lambda : w_allowed[random.random_integers(0,2)]
        relt = lambda x: x[random.random_integers(0,len(x) -1)]
        gelt = lambda x: g.append(relt(x)) or ga.append(x)

        print randwt()

        ga.append(n_hc_allowed)
        na = relt(n_hc_allowed) #only one element to select...
        g.append(na)
        for i in range(na):
            gelt(hc_allowed)
            gelt(hcc_allowed)
            gelt(w_allowed)
        
        ga.append(n_hgs_allowed)
        nhgs = relt(n_hgs_allowed)
        g.append(nhgs)
        for i in range(nhgs):
            gelt(tg_allowed)
            gelt(hgs_allowed)
            gelt(w_allowed)
        
        ga.append(n_hgu_allowed)
        nhgu = relt(n_hgu_allowed)
        g.append(nhgu)

        #NOTE!
        #HGU SHOULD BE ENCODED AS
        #A VECTOR WHERE DIFFERENT
        #GENES HAVE DIFFERENT U's 
        #ALLOWED... FOR NOW, JUST
        #HAVE HGU_ALLOWED = MAX_U
        for i in range(nhgu):
            gelt(tg_allowed)
            gelt(hgu_allowed)
            gelt(w_allowed)

            
        ga.append(n_ogc_allowed)
        nogc = relt(n_ogc_allowed)
        g.append(nogc)
        for i in range(nogc):
            gelt(ogc_allowed)
            gelt(ogcc_allowed)
            gelt(w_allowed)
        
        for i in range(self.ntg):
            gelt(w_allowed)
        
        print g
        return g, ga
        
  def make_cxns_from_genome(self,genome):
    #reverse the genome to read through by popping
        genome = genome[::-1]

        hc_cxns = [[]] * self.nhc
        hc_ws = [[]] * self.nhc
        hgs_cxns = [[]] * self.ntg
        hgs_ws = [[]] * self.ntg
        hgu_cxns = [[]] * self.ntg
        hgu_ws = [[]] * self.ntg
        ogg_ws = [[]] * self.ntg
        ogc_cxns = [[]] * self.ntg
        ogc_ws = [[]] * self.ntg

        #connections between cluster hidden nodes and shared tfs
        #n = ns  * 2?
        #for the moment, numbers are not being allowed to vary...
        #that would create big problems, I think.
        n_hc_cxns = genome.pop()
        for i in range(n_hc_cxns):
            idx, cxn, w = [genome.pop() for j in range(3)]
            if cxn in hc_cxns[idx]: continue
            hc_cxns[idx].append(cxn)
            hc_ws[idx].append(w)
   
        #connections between tg hidden nodes and shared tfs
        #n = 0?
        n_hgs_cxns = genome.pop()
        for i in range(n_hgs_cxns):
            idx, cxn, w = [genome.pop() for j in range(3)]
            if cxn in hgs_cxns[idx]:continue
            hgs_cxns[idx].append(cxn)
            hgs_ws[idx].append(w)
            
        
        #connections between tg hidden nodes and unshared tfs
        #n ~ 3?
        n_hgu_cxns = genome.pop()
        for i in range(n_hgu_cxns):
            idx, cxn, w = [genome.pop() for j in range(3)]
            if cxn in hgu_cxns[idx]:continue
            hgu_cxns[idx].append(cxn)
            hgu_ws[idx].append(w)

        
        #connection between tg outputs and cluster hidden nodes
        #n = nc * ng?
        n_ogc_cxns = genome.pop()
        for i in range(n_ogc_cxns):
            idx, cxn, w = [genome.pop() for j in range(3)]
            if cxn in ogc_cxns[idx]:continue
            ogc_cxns[idx].append(cxn)
            ogc_ws[idx].append(w)

        #NOTE THAT THE NUMBER OF OGG CONNECTIONS = NG
        #(not read from genome)
        n_ogg_cxns = self.ntg
        for i in range(n_ogg_cxns):
            ogg_ws[i] = genome.pop()

        assert(genome == [])
    
        self.mynn.set_connections(hc_cxns,
                                 hc_ws,
                                 hgs_cxns,
                                 hgs_ws,
                                 hgu_cxns,
                                 hgu_ws,
                                 ogg_ws,
                                 ogc_cxns,
                                 ogc_ws)
        
  def init_net(self):
            
        ng = self.ntg
        nc = self.nhc
        nx = self.ntf_s + sum(self.ntf_u)
        nh = ng + nc

        in_mod_sizes = [self.ntf_s] + self.ntf_u
        h_mod_sizes =  [self.nhc] + [1 for i in range(self.ntg)]
        out_mod_sizes = [1 for i in range(self.ntg)]
        self.mynn = MyNN(in_mod_sizes, h_mod_sizes, out_mod_sizes)
    
  def set_data(self,xdat, ydat):
        self.mynn.set_data(xdat,ydat)
        
  def set_trainer(self):
        self.mynn.set_trainer()
  def train(self):
        self.mynn.train()

    




class Net():
    def __init__(self, nin, nh, nout ):

        from pybrain.structure import FeedForwardNetwork
        net = FeedForwardNetwork()
        

        self.nin = nin
        self.nh = nh
        self.nout = nout

        self.net = net

    def set_data(self,xdat,ydat):
        if shape(xdat)[0] != shape(ydat)[0]:
            raise Exception('dimension mismatch b/w x, y')
        nt = len(xdat)
        ny = shape(ydat)[1]
        nx = shape(xdat)[1]

        from pybrain.datasets import SupervisedDataSet
        ds = SupervisedDataSet(nx, ny)
        for i in range(nt):
            ds.addSample(xdat[i], ydat[i])
        self.data = ds

    def set_connections(self,
                        hc_cxns,
                        hc_ws,
                        hgs_cxns,
                        hgs_ws,
                        hgu_cxns,
                        hgu_ws,
                        ogg_ws,
                        ogc_cxns,
                        ogc_ws):


        from pybrain.structure import FullConnection
        net = self.net

        ins = []
        hids = []
        outs = []

        for i in range(len(self.nin)):
            ins.append(LinearLayer(self.nin[i],name = 'i'+str(i)))
        for l in ins:
            self.net.addInputModule(l)

        for i in range( len(self.nh)):
            hids.append(SigmoidLayer(self.nh[i],name = 'h'+str(i)))

            print 'nh: ',self.nh[i]
        for l in hids:
            print 'len ',l.dim
            self.net.addModule(l)

        for i in range( len(self.nout)):
            outs.append(LinearLayer(self.nout[i], name = 'o'+str(i)))
        for l in outs:
            self.net.addOutputModule(l)

        #connect up the hidden clusters.
            
        name_wts = {}
        for i in range(len(hc_cxns)):
            for j in range(len(hc_cxns[i])):
                cstart = hc_cxns[i][j]
                cw = hc_ws[i][j]
                inmod = ins[0]
                instart = cstart
                infinish = instart+1
                outmod = hids[0]
                outstart = i
                outfinish = outstart + 1

                name = 'hc'+str(i)+str(j)
                name_wts[name] = cw 

                cxn = FullConnection(inmod, outmod, 
                                     inSliceFrom = instart,inSliceTo = infinish,
                                     outSliceFrom = outstart, outSliceTo = outfinish,
                                     name = name)
                net.addConnection(cxn)
            
 

            
        for i in range(len(hgs_cxns)):
            for j in range(len(hgs_cxns[i])):
                cstart = hgs_cxns[i][j]
                cw = hgs_ws[i][j]
                inmod = ins[0]
                instart = cstart
                infinish = instart +1
                outmod = hids[i+1]
                
                #only one elements in each module...
                outstart = 0
                outfinish = outstart +1

                name = 'hgs'+str(i)+str(j)
                name_wts[name] = cw 


                cxn = FullConnection(inmod, outmod, 
                                     inSliceFrom = instart,inSliceTo = infinish,
                                     outSliceFrom = outstart, outSliceTo = outfinish,
                                     name = name)
                #cxn.params[0] = cw
                net.addConnection(cxn)

  
             
        for i in range(len(hgu_cxns)):
            for j in range(len(hgu_cxns[i])):
                cstart = hgu_cxns[i][j]
                cw = hgu_ws[i][j]
                inmod = ins[i+1]
                
                #NOTE, I HAVE TO MOD THE HGUs BECAUSE 
                #HGU CAN RANGE FROM 0...MAX_HGU in the genome
                instart = mod(cstart,inmod.dim)
                infinish = instart +1
                outmod = hids[i+1]
                

                #only one elements in each module...
                outstart = 0
                outfinish = outstart + 1 
                
                name = 'hgu'+str(i)+str(j)
                name_wts[name] = cw 

                cxn = FullConnection(inmod, outmod, 
                                     inSliceFrom = instart,inSliceTo = infinish,
                                     outSliceFrom = outstart, outSliceTo = outfinish,
                                     name = name)
                #cxn.params[0] = cw

                print name
                print instart,infinish,outstart,outfinish
                print inmod.dim, outmod.dim
                print inmod.name, outmod.name
                net.addConnection(cxn)  

             
        for i in range(len(ogc_cxns)):
            for j in range(len(ogc_cxns[i])):
                cstart = ogc_cxns[i][j]
                cw = ogc_ws[i][j]
                inmod = hids[0]
                instart = cstart
                infinish = instart + 1
                outmod = outs[i]
                
                #only one elements in each module...
                outstart = 0
                outfinish = outstart +1
                                
                name = 'ogc'+str(i)+str(j)
                name_wts[name] = cw 

                cxn = FullConnection(inmod, outmod, 
                                     inSliceFrom = instart,inSliceTo = infinish,
                                     outSliceFrom = outstart, outSliceTo = outfinish,
                                     name = name)


                net.addConnection(cxn)
  
        for i in range(len(ogg_ws)):
            inmod = hids[(i + 1)]
            outmod = outs[i]
            cxn = FullConnection(inmod,outmod,name =name)
            #cxn.params[0] = cw
                                            
            name = 'ogg'+str(i)
            name_wts[name] = cw 

            net.addConnection(cxn)
            

        self.net.sortModules()

        for v in net.connections.items():
            for elt in v[1]:
                name = elt.name
                elt.params[0] = name_wts[name]
        return        



    def set_trainer(self):        
        trainer = BackpropTrainer(self.net, self.data)
        self.trainer = trainer

    def train(self):
        if not self.trainer or not self.data: 
            raise Exception()
        for t in range(10):
            t = self.trainer.train()
     	


class Pop():
  pass

#DEMO CLASS FROM PYBRAINZ DOC
class MyPB():

    def run(self):
        genome = G1DList.G1DList(20)
        genome.evaluator.set(eval_func)
        ga = GSimpleGA.GSimpleGA(genome)
        ga.evolve(freq_stats=10)
        print ga.bestIndividual()


