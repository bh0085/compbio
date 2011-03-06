#!/usr/bin/env python
import compbio.projects.cbdb.ncbi_tax_gbjoin as gbj
if __name__ == '__main__':
  print 'filling db'
  gbj.fill_db()
  exit(0)
