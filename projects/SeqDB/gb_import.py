from Bio import SeqIO
import os
import compbio.config as config
from SeqDB import *
from sqlalchemy.orm.exc import NoResultFound
from numpy import *

def add_archaea():
  gb_paths =[] 
  for root, dirs, files in os.walk(os.path.join(config.root, 'sequences/16s')):
    for f in files:
      if '.gbk' in f:
        gb_paths.append(os.path.join(root, f))

  gb_paths.remove(gb_paths[0])
  skip = True
  for p in gb_paths:
    print 'importin path: '
    print '   ' + p
    if  '247939seqs' in p: 
      skip = False
    if skip or '247939seqs' in p:
      continue

    import_alignment( p)

def import_alignment( filename, alignment_name = None):
  f = open(filename)
  if not alignment_name:
    alignment_name = os.path.splitext(os.path.basename(filename))[0]
  alignment = tryAddUnique(Alignment(alignment_name))
  recs =  SeqIO.parse(f, 'genbank')
  count = 0
  

  for rec in recs:
    count += 1
    if mod(count, 1000) == 0: print count
    #if count <    169000: continue 

    oname = rec.annotations['organism']
    tlist = rec.annotations['taxonomy']
    gbid   = rec.annotations['comment'][9:]
    seq   = rec.seq
    organism = rec.annotations['organism']
    seq_source = rec.annotations['source']
    seq_description = rec.description
    seq_name = rec.name
    seq_entry = tryAddUnique(Sequence(seq.__str__(), seq_name, seq_source, seq_description),
                             name = seq_name)
    alignment_join = tryAddUnique(Alignment_SequenceJoin(alignment.id,
                                                         seq_entry.id),
                                  sequence = seq_entry.id)
    gbid_entry = tryAddUnique(GBID(gbid),name = gbid)
    gbid_join = tryAddUnique(GBID_SequenceJoin(gbid_entry.id,
                                               seq_entry.id),
                             sequence = seq_entry.id)
    gb_entry = tryAddUnique(GBID(gbid), name = gbid)
    tax_units = [Domain, Phylum, Class, Order, Family, Genus]
    tax_joins = [Domain_GBIDJoin, Phylum_GBIDJoin, Class_GBIDJoin, \
                   Order_GBIDJoin, Family_GBIDJoin, Genus_GBIDJoin \
                   ]

    if len(tlist) == 7:
      tlist = tlist[1:7]
    if 'Bacteria' in tlist[1]:
      if len(tlist) >4:
        tlist = tlist[0:4]+tlist[5:]
      if len(tlist)>5:
        tlist = tlist[0:5]+tlist[6:]
    tlist = tlist[1:]
    for i in range(len(tlist)) :
      u = tax_units[i]
      tname = tlist[i].replace('"','')
      tax_entry = tryAddUnique(u(tname), name = tname)
      tax_gbid_join = tryAddUnique(tax_joins[i](tax_entry.id,gb_entry.id) ,
                                   gbid = gb_entry.id)
      
                                   
    spec = Specie(organism)
    spec_entry =tryAddUnique(spec, name = organism)
    spec_gbid_join = tryAddUnique(Specie_GBIDJoin(spec_entry.id, 
                                                  gb_entry.id),
                                  gbid = gb_entry.id)

    if mod(count, 1000) == 0:
      Session.commit()
      print 'commiting, ',count, tax_entry.name
    

  print 'commiting, done with file {0}'.format(filename)
  Session.commit()
