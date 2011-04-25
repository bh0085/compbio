#!/usr/bin/env python
'''
bsruns.py

Scripts designed to be called from the commandline with a
run ID (as opposed to a bsub jobID) as the 

'''

from compbio.utils.bsub_utils import *

def bout(run_id):
    '''
Returns the output of the program having a given runid.

usage:      bout(run_id)
commandline bsruns.py bout run_id

'''
    return load_data(run_id, 'output') 

def bstatus(run_id):
    '''
Get the status of a currently running job.
By convention, statuses will be saved in json
format so that bstatus is compatible with printing
to the stdout and communication over ssh.

usage:       bstatus(run_id)
commandline: bsruns.py bstatus run_id
'''
    return load_data(run_id, 'status')


if __name__ == '__main__':
    assert len(sys.argv) > 2
    if sys.argv[1] in ['bout', 'bstatus' ]:
        #Dump data to json and write to stdout
        run_id = sys.argv[2]
        sys.stdout.write(sjson.dumps(globals()[sys.argv[1]](run_id)))
        exit(0)
    else:
        raise Exception()
