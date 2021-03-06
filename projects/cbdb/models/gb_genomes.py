from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import relation
import compbio.projects.cbdb
from compbio import config
import os
import numpy as np
from Bio import SeqIO
from itertools import ifilter

def get_tables():
    #by convention, store alignment files in dataPath/alignments.
    return [
            dict(name = 'Genome',
                 attrs={'id':Column(Integer, primary_key = True),
                        'name':Column(String),
                        'seq_url':Column(String),
                        'source_taxon':Column(Integer,index = True),
                        'source_organism':Column(String),
                        'gb_accession':Column(String),
                        'annotations':Column(String),
                        }),
            dict(name = 'Feature',
                 attrs={'id':Column(Integer, primary_key = True),
                        'start':Column(Integer),
                        'start_ext':Column(Integer),
                        'end':Column(Integer),
                        'end_ext':Column(Integer),
                        'type':Column(String),
                        'strand':Column(Integer),
                        'genomeid':Column(Integer, ForeignKey('genome.id')),
                        'genomeobj':relation('Genome', backref = 'features')
                        }),
            dict(name = 'SubFeature',
                 attrs={'id':Column(Integer, primary_key = True),
                        'start':Column(Integer),
                        'start_ext':Column(Integer),
                        'end':Column(Integer),
                        'end_ext':Column(Integer),
                        'type':Column(String),
                        'strand':Column(Integer),
                        'featureid':Column(Integer, ForeignKey('feature.id')),
                        'featureobj':relation('Feature', backref = 'subfeatures')}),
            dict(name = 'Qualifier',
                 attrs = {'id':Column(Integer, primary_key = True),
                          'key':Column(String),
                          'value':Column(String),
                          'featureid':Column(Integer, ForeignKey('feature.id')),
                          'featureobj':relation('Feature',backref = 'qualifiers'),
                          'subfeatureid':Column(Integer, ForeignKey('subfeature.id')),
                          'subfeatureobj':relation('SubFeature',backref='qualifiers')}
                 )
            ]
def fill_db( name = 'bacterial_genomes', reset = False,
              postgres = False, host = 'broad'):
    dbi = cbdb.getName(
                       name,
                       postgres = postgres,
                       tables = get_tables(),
                       reset = np.mod(reset, 2), 
                       host = host)


    paths = []
    for r,ds, fs in os.walk('/Volumes/ganymede/all.gbk/'):
      for f in fs:
        if 'gbk' in f: paths.append(os.path.join(r, f))
        count = 0 
    

    for p in paths:
      
      if count < 1668:
        count += 1
        continue
      count += 1
      fopen = open(p)
      for rec in SeqIO.parse(fopen, 'genbank'):
        f0 = rec.features[0]
        if f0.type == 'source':
          source_taxon = f0.qualifiers['db_xref'][0][6:]
          source_organism=f0.qualifiers['organism'][0]
        else:
          source_taxon = None
          source_organism = None
          
        fa_seqpath = 'genomes/'+rec.id+'.fa'
        fa_sequrl = config.dataURL(fa_seqpath)
        fa_seqfile = config.dataPath(fa_sequrl)
        fopen = open(fa_seqfile,'w')
        SeqIO.write(rec,fopen, 'fasta')
        fopen.close()

        adds = []
        genome = dbi.Genome(name = rec.name, 
                           seq_url =fa_sequrl,
                           source_taxon = source_taxon,
                           source_organism = source_organism,
                           gb_accession = rec.id,
                           annotations = rec.annotations.__str__())

        #adds.append(genome)
        print 'adding genome ' + source_organism
        dbi.Session.add(genome)
        print 'commiting update ' 
        dbi.Session.commit()
        print 'genome added! '
        for f in rec.features:
          feature = dbi.Feature(type = f.type,
                                start = f.location.start.position,
                                start_ext = f.location.start.extension,
                                end = f.location.end.position,
                                end_ext = f.location.end.extension,
                                strand = f.strand,
                                genomeobj = genome)
          #print 'adding feature ' + f.type
          #dbi.Session.add(feature)
          adds.append(feature)
          for k,v in f.qualifiers.iteritems():
            q = dbi.Qualifier(key = k,
                                      value = v.__str__(),
                                      featureobj = feature)
            #dbi.Session.add(q)
            adds.append(q)
          for sf in f.sub_features:
            sub = dbi.SubFeature(type = sf.type,
                                 start = sf.location.start.position,
                                 start_ext = sf.location.start.extension,
                                 end =sf.location.end.position,
                                 end_ext = sf.location.end.extension,
                                 strand = sf.strand,
                                 featureobj = feature)
            adds.append(sub)
            #dbi.Session.add(sub)
            for k,v in sf.qualifiers.iteritems():
              q = dbi.Qualifier(key = k,
                                value = v.__str__(),
                                subfeatureobj = sf)
              #Session.add(q)
              adds.append(q)
                                
        dbi.Session.add_all(adds)



        if np.mod(count, 2) == 0:
          print count
#print count, p , seq.source_organism
          print 'committing update'
          dbi.Session.commit()
          print 'update commited!'
      dbi.Session.commit()
    
