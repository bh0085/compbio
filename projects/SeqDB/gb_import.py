from Bio import SeqIO
import os
import compbio.config as config
from sqlalchemy.orm.exc import NoResultFound
from numpy import *
import compbio.projects.SeqDB as sqd

def add_archaea(db):
  gb_paths =[] 
  for root, dirs, files in os.walk(os.path.join(config.root, 'sequences/16s')):
    for f in files:
      if '.gbk' in f:
        gb_paths.append(os.path.join(root, f))

  gb_paths.remove(gb_paths[0])
  skip = True
  for p in gb_paths:
    import_alignment( p, db)

def import_alignment( filename, db,  alignment_name = None):
  f = open(filename)
  if not alignment_name:
    alignment_name = os.path.splitext(os.path.basename(filename))[0]
  alignment = db.tryAddUnique(db.Alignment(alignment_name))
  gbfile = db.tryAddUnique(db.Gbfile(filename))
  recs =  SeqIO.parse(f, 'genbank')
  count = 0
  

  for rec in recs:
    count += 1
    oname = rec.annotations['organism']
    tlist = rec.annotations['taxonomy']
    gbid   = rec.annotations['comment'][9:]
    seq   = rec.seq
    organism = rec.annotations['organism']
    seq_source = rec.annotations['source']
    seq_description = rec.description
    seq_name = rec.name
    seq_entry = db.tryAddUnique(db.Sequence(seq.__str__(), seq_name, seq_source, seq_description),
                             name = seq_name)
    alignment_join = db.tryAddUnique(db.Alignment_SequenceJoin(alignment.id,
                                                         seq_entry.id),
                                  sequence = seq_entry.id)
    gbid_entry = db.tryAddUnique(db.GBID(gbid, f.tell()),name = gbid)
    gbid_join = db.tryAddUnique(db.GBID_SequenceJoin(gbid_entry.id,
                                                     seq_entry.id),
                            sequence = seq_entry.id)
    file_join = db.tryAddUnique(db.Gbfile_GBIDJoin(gbfile.id,gbid_entry.id),
                                gbid = gbid_entry.id)
    tax_units = [db.Domain, db.Phylum, db.Class, db.Order, db.Family, db.Genus]
    tax_joins = [db.Domain_GBIDJoin, db.Phylum_GBIDJoin, db.Class_GBIDJoin, \
                   db.Order_GBIDJoin, db.Family_GBIDJoin, db.Genus_GBIDJoin \
                   ]

    if 'Bacteria' in tlist[1]:
      if len(tlist) >4:
        tlist = tlist[0:4]+tlist[5:]
      if len(tlist)>5:
        tlist = tlist[0:5]+tlist[6:]
    tlist = tlist[1:]
    for i in range(len(tlist)) :
      u = tax_units[i]
      tname = tlist[i].replace('"','')
      tax_entry = db.tryAddUnique(u(tname), name = tname)
      tax_gbid_join = db.tryAddUnique(tax_joins[i](tax_entry.id,gbid_entry.id) ,
                                   gbid = gbid_entry.id)
      
                                   
    spec = db.Specie(organism)
    spec_entry =db.tryAddUnique(spec, name = organism)
    spec_gbid_join = db.tryAddUnique(db.Specie_GBIDJoin(spec_entry.id, 
                                                  gbid_entry.id),
                                  gbid = gbid_entry.id)

    if mod(count, 1000) == 0:
      db.Session.commit()
      print 'commiting, ',count, tax_entry.name
    if mod(count, 5000) == 1001:
      print 'domains so far: '
      print map(lambda x: x.name, db.Session.query(db.Domain).all())
      print 'phyla so far: '
      print map(lambda x: x.name, db.Session.query(db.Phylum).all())
    

  print 'commiting, done with file {0}'.format(filename)
  db.Session.commit()
