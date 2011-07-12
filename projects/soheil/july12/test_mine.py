#!/usr/bin/env python
import cb.config as cfg
import os
sdp = cfg.dataPath('soheil')
matfiles = [f for f in os.listdir(sdp) if '.mat' in f]

for f in matfiles: 
    infile, resfile, totfile = \
        f, \
        os.path.join(*( os.path.split(f)[:-1]+('mat_out',)+os.path.split(f)[-1:])),\
        os.path.join(*( os.path.split(f)[:-1]+('mat_all_out',)+os.path.split(f)[-1]))
    print infile, resfile, totfile
    

exit
