from Bio import SeqIO
import os, re
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
  
