#!/usr/bin/env python
import compbio.config as cfg
import sys

def usage():
	print '''
Usage: 
  pydatapath.py path volume
'''

if __name__ == '__main__':
	path = sys.argv[1] if len(sys.argv) > 1 else ''
	volume = sys.argv[2] if len(sys.argv) > 2 else ''
	host = ''
	
	if path == 'usage':
		usage()
		exit(0)
	

	localpath = sys.argv[1] if len(sys.argv) > 1 else ''
	path = cfg.dataPath(cfg.dataURL(path, volume_name = volume))
	sys.stdout.write(path)
	exit(0)	
