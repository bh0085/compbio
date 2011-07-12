#!/usr/bin/env python
import cb.config as cfg
sdp = cfg.dataPath('soheil')
matfiles = [f for f in os.listdir(sdp) if '.mat' in f]

for f in matfiles: 
    print f

exit
