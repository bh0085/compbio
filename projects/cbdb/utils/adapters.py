import os
import numpy as np

def fileMap2DB(dbname,
               table_files,
               table_names,
               table_defs,
               column_maps,
               col_sep,
               rec_sep,
               starts = None,
               finishes = None,
               reset = False,
               ncommit = 10000
               ):
  rec_iterfun = lambda x: x.xreadlines()
  dbi = cbdb.getName(dbname, 
                     tables = table_defs,
                     reset =reset)
  
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
    if start != None:
      fopen.seek(start)
    for l in rec_iterfun(fopen):
      if (finish != None):
        if fopen.tell() >= finish:
          break
      #read ahead until the record seperator is encountered.
      count += 1
      l0 += l
      if rec_sep in l0:
        l = l0
        l0 = ''
      else: continue
      
      cols = colfun(l)
      cls = mapped_class(**dict(map(lambda (x,y):(x,cols[y]),
                                    mapped_columns.iteritems())))
      dbi.Session.merge(cls)
      if np.mod(count, ncommit) == 0:
        dbi.Session.commit()
        print k,v, count, cols
    dbi.Session.commit()
  return

def fileMap2DB2(dbname,
               table_files,
               table_names,
               table_defs,
               column_maps,
               col_sep,
               rec_sep,
               read_inc,
               reset = False,
               ncommit = 10000
               ):
  rec_iterfun = lambda x: x.xreadlines()
  dbi = cbdb.getName(dbname, 
                     tables = table_defs,
                     reset =reset)

  
  joined = [(table_names[i],table_files[i]) \
                               for i in range(len(table_names))]
  fill_tables = dict(joined)
  colfun = lambda x: unicode(x, errors = 'replace').replace(rec_sep, '').split(col_sep)

  count = 0
  for k, v, in fill_tables.iteritems():
    fopen = open(v)
    maxsize = os.path.getsize(v)
    mapped_class = dbi.__dict__[k]
    mapped_columns = column_maps[k]

    while fopen.tell() + 1 < maxsize:
      lines =  fopen.readlines(read_inc)
      cols = map(lambda x: colfun(x), lines)
      cls = [mapped_class(**dict(map(lambda (x,y):(x,c[y]),
                                     mapped_columns.iteritems())))
             for c in cols]

      dbi.Session.add_all(cls)
      print 'commiting ' + str(len(lines))
      dbi.Session.commit()
      
      print 'done committing, through with %5.3f ' % ((float(fopen.tell()))/maxsize)
  return

import compbio.projects.cbdb
