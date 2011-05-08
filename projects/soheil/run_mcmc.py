#!/usr/bin/env python
import sys
import compbio.utils.bsub_utils as bsu

def launch_many(run_id):
    '''
Generate script paramaters and launch a bunch of bsub jobs.

Designed to be run on the cluster via an interactive shell.
'''
    print 'Launching all jobs!'


def run_single(run_id):
    '''
Given an input dictionary containing a single paramater set
run mcmc in matlab using bs_macros.run_matlab.
'''



if __name__ == '__main__':
    run_id = sys.argv[2]
    run_func = globals()[sys.argv[1]]
    output_dict = run_func(run_id)
    if output_dict == None:
        output_dict = {'blank':'Nothing output in call to {0}'.\
                           format(sys.argv[1])}
    bsu.save_data( output_dict, run_id, 'output')
