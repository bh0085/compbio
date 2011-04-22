#!/usr/bin/env python
import compbio.config as cfg
import sys

if __name__ == '__main__':
	localpath = sys.argv[1] if len(sys.argv) > 1 else ''
	path = cfg.dataPath(localpath)
	sys.stdout.write(path)
	exit(0)	
