from compbio.projects import cbdb
import os
import numpy as np

def fileMap2DB(dbname,
               table_files,
               table_names,
               table_defs,
               column_maps,
               col_sep,
               rec_sep,
               reset = False
               ):
  rec_iterfun = lambda x: x.xreadlines()
  dbi = cbdb.getName(dbname, 
                     tables = table_defs,
                     reset = np.mod(reset, 2))
  
  joined = [(table_names[i],table_files[i]) \
                               for i in range(len(table_names))]
  fill_tables = dict(joined)
  colfun = lambda x: unicode(x, errors = 'replace').replace(rec_sep, '').split(col_sep)

  count = 0
  for k, v, in fill_tables.iteritems():
    fopen = open(v)
    mapped_class = dbi.__dict__[k]
    mapped_columns = column_maps[k]
    
    l0 = ''
    for l in rec_iterfun(fopen):
      #read ahead until the record seperator is encountered.
      count == 1
      l0 += l
      if rec_sep in l0:
        l = l0
        l0 = ''
      else: continue
      
      cols = colfun(l)
      cls = mapped_class(**dict(map(lambda (x,y):(x,cols[y]),
                                    mapped_columns.iteritems())))
      dbi.Session.merge(cls)
      if np.mod(count, 1000) == 0:
        dbi.Session.commit()
        print k,v, count, cols
    dbi.Session.commit()
  return
