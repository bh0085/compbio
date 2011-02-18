from Bio import SeqIO
import os, re, itertools as it
from compbio import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import compbio.utils.hilbert as h, utils.plots as plots, utils.colors as mycolors
from compbio.utils.path_mgr import *
def mycGB():
  dobr = 0
  for root, dirs, files in os.walk('/data/genomes/Bacteria'):
    for f in files:
      if re.search(re.compile('gbk$'), f):
        g = os.path.join(root, f )
        gb_record = SeqIO.read(open(g,"r"), "genbank")
        print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
        print repr(gb_record.seq)
        
        for f in gb_record.features[0:10]:
          print f
          #break
          dobr = 1
    if dobr:
      break

def align():
  pass


def findseq():
  gpath = config.root
  mycgbk = os.path.join(gpath, 'genomes/Bacteria/Mycoplasma_mycoides_SC_PG1_uid58031/NC_005364.gbk')
  mycfa = os.path.join(gpath, 'genomes/Bacteria/Mycoplasma_mycoides_SC_PG1_uid58031/NC_005364.fna')
  grec = SeqIO.parse(mycgbk, 'genbank')
  frec = SeqIO.parse(mycfa, 'fasta')
  fa = list(frec)[0]

  grecf0 = list(grec)[0].features
  grecf = it.ifilter(lambda x: 'gene' in x.qualifiers.keys(), grecf0)
  grecf = sorted(list(grecf), key = lambda x:x.qualifiers['gene'])

  geneproducts = {}

  
  for key, grp in it.groupby(grecf, lambda x : x.qualifiers['gene']):
    g = list(grp)
    print key, len(g)

    for k in key:
      geneproducts[k] = g

  lens =  map(lambda x: len(x), geneproducts.values())
    
  loc = lambda x: [x.location.start.position, x.location.end.position]
  expand = lambda x, y: [np.min(x[0],y[0]), np.max(x[1],y[1])]
                    
  nprods = len(geneproducts.keys())
  import compbio.utils.colors
  colors = utils.colors.getct(nprods)
  cshow = range(10)
  clam = lambda x: x in cshow and colors[x]  or [1,0,0]
  
  
  f = plt.figure(0)
  f.clear()
  ax = f.add_subplot(111)
  rs, cs = [],[]
  ct = 0


  #utils.plots.spacefill(ax, rs, cs)
  import compbio.utils.hilbert as h

  speed_hash = h.gInterp(grecf0)

  #levs = h.inverseN(len(rs))
  bert = h.Hilbert(5, angle = 90)
  all_verts = bert.vertices(0,1)

  
  ax.plot(all_verts[:,0],
          all_verts[:,1],
          color = 'black',
          alpha = .1)
  
  src = list(it.ifilter(lambda x: x.type =='source',grecf0))[0]
  src_pos = [src.location.start.position,
             src.location.end.position]
  genome_l = src_pos[1] - src_pos[0]
  
  for i in range(len(geneproducts.keys())):
    g=geneproducts.values()[i]
    for elt in g:
      if elt.type == 'CDS':
        r = array([elt.location.start.position , elt.location.end.position],float)
        fracs = [speed_hash[r[0]],speed_hash[r[1]]]
        verts = bert.vertices(fracs[0],fracs[1],snap = False)
        
        ax.plot(verts[:,0],verts[:,1], 
             color = clam(i),
             linewidth = 8)

  raise Exception()



def testBlast():
  from Bio.Blast.Applications import NcbiblastxCommandline
  #help(BlastallCommandline)
  dpath = os.path.join(config.root, 'genomes/myco.fna')
  qpath = os.path.join(config.root, 'genomes/test.txt')
  opath = os.path.join(config.root, 'genomes/out.xml')

  #blastall -p blastn -d myco.fna -i test.txt -o test.r

  blastx_cline = NcbiblastxCommandline( query=qpath, db=dpath, 
                                     out = opath)
  print blastx_cline
  stdout, stderr = blastx_cline()
  
  raise Exception()
  stdout, stderr = blastx_cline()

  
  raise Exception()

def getGBKS():
  gpath = os.path.join(config.root, 'genomes/Bacteria')
  
  gbks = []
  for root, folders, files in os.walk(gpath):
    for f in files:
      if re.search('.*gbk', f):
        gbks.append(os.path.join(root, f))
    if len(gbks) >= 1:
      break

  nx = int(ceil(sqrt(len(gbks))))
  ny = int(ceil(sqrt(len(gbks))))
  
  f = plt.figure(1)

  f.clear()
  for i in range(len(gbks)):
    xv,yv = divmod(i,nx)
    ax = f.add_subplot('{0}{1}{2}'.format(nx,ny, i), xlim = [-.2,1.2], ylim = [-.2,1.2], frameon=True)
    
    maptRNAs(gbks[i], ax)


tcols = {}
tct =  mycolors.getct(64)
lvls = 7
functype = 'blank'

def maptRNAs(gbkfile, ax = None):

  grec = SeqIO.parse(gbkfile, 'genbank')

  grec = list(grec)[0].features
  blank_speedfun = lambda x :\
      x.type == 'tRNA' and 2 \
      or x.type == 'gene' and 2 \
      or x.type == 'CDS' and 2 \
      or 0

  speedfun = lambda x :\
      x.type == 'tRNA' and 5000 \
      or x.type == 'gene' and 1 \
      or x.type == 'CDS' and 5 \
      or 0

  global funtype
  if functype == 'speed':
    speed_hash = h.gInterp(grec,speedfun)
  else:
    speed_hash = h.gInterp(grec,blank_speedfun)

  global lvls
  bert = h.Hilbert(lvls, angle = 90)
  
  color_lam = lambda x: x.type =='tRNA' and 'red' \
      or x.type =='CDS' and 'blue' \
      or 'none'
  
  z_lam = lambda x: x.type =='tRNA' and '1' \
      or 0
  
  width_lam = lambda x: x.type =='tRNA' and 8 \
      or 16
  
  alpha_lam = lambda x: x.type =='tRNA' and 1 \
      or .5


  if not ax:
    f = plt.figure(0)
    f.clear()
    ax = f.add_axes([0,0,1,1], xlim = [-.2,1.2],
                    ylim = [-.2,1.2])

  from matplotlib.collections import LineCollection as lc
  

  v0 = bert.vertices(0,1)
  ax.plot(v0[:,0], v0[:,1], color = 'black', alpha = .2)

  global tcols
  global tct
  for r in grec:
    if r.type != 'tRNA':
      continue
    p = r.qualifiers['product'][0]
    aa = re.search(re.compile('-(.*)'), p).group(1)
    if not aa in tcols.keys():
      tcols[aa] = tct.pop()

  print tcols['Phe']

  for r in grec:
    if not (r.type == 'gene' or r.type =='CDS' or r.type == 'tRNA'):
      continue

    color = 'none'
    if r.type == 'tRNA':
      p = r.qualifiers['product'][0]
      aa = re.search(re.compile('-(.*)'), p).group(1)
      color = tcols[aa] 

    if r.type == 'gene' or r.type == 'CDS':
      rec = r
      string = rec.__str__()

      for k, v in tcols.iteritems():
        if k in string and (not k in ['Val','Pro']) :
          color = v          
            
      if 'ynthetase' in string:
        color = 'black'
      
        
    if color == 'none':
      continue



    #color = color_lam(r)
    alpha = alpha_lam(r)
    width = width_lam(r)
    z = z_lam(r)

    l0 = r.location.start.position
    l1 = r.location.end.position
    fracs = speed_hash[l0],speed_hash[l1]
    verts = bert.vertices(fracs[0],fracs[1], snap = False)

    ax.plot(verts[:,0], verts[:,1],
            color = color,
            linewidth = width,
            alpha = alpha,
            zorder = z)
    


def compileDB():
  gdir = os.path.join(config.root, 'genomes/Bacteria')
  frecs = []
  grecs = []
  for root, dirs, files in os.walk(gdir):
    for f in files:
      if re.search(re.compile('fna$'), f):
        fr = SeqIO.read(open(os.path.join(root, f)),'fasta')
        frecs.append(fr)
        
      if re.search(re.compile('gbk$'), f):
        fr = SeqIO.read(open(os.path.join(root, f)),'genbank')
        grecs.append(fr)

        
        
  fout = os.path.join(config.root,'genomes/myco.fna')
  gout = os.path.join(config.root,'genomes/myco.gbk')
  SeqIO.write(frecs, fout, 'fasta')
  SeqIO.write(grecs, gout, 'genbank')

def bpTut():  
  fa = open(os.path.join(config.root, 'sequences/ls_orchid.fasta'))
  gb = open(os.path.join(config.root, 'sequences/ls_orchid.gbk'))
  for record in SeqIO.parse(fa,'fasta'):
    print record
  for g in SeqIO.parse(gb, 'genbank'):
    for f in g.features:
      print f
    
