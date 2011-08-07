#!/usr/bin/env bash

'''
Take the ideas from day1 and run them in a grid
search fashion on the broad. Now that we are using
the cluster, can do an l^2 analysis.

For 1d analysis:
run0(method = '1d')


'''

import day0
from numpy import *

def run0(method = '1d', win_len = 150, win_ofs = 25):
    '''
Build a list of jobs on the Broad within which 
to compute zscores for every single locus in the
PVT1 region.

'''
    ali_nums = day0.fetch_num_ali()
    runs_per_rng = 1e5 / 25 
    
    p0 = 0; p = []
    while p0 + win_len < shape(ali_nums)[1]:
        b0 = p0
        p0 += win_ofs * runs_per_rng
        b1 = p0
            
        
        p.append(dict(win_len = win_len,
                      win_ofs = win_ofs,
                      baserng = (b0, b1),
                      run_id = 'run1d_single_{0}_{1}'))

    runid = 'run1d_{0}_{1}'.format(win_len,win_ofs)
    ll = run1_launcher(p, runid)
    ll.launch()

    return ll
         
                 

#Local launchpoint for bsubs
def run1d_launcher(input_dicts, run_id, host = 'tin' ):
    '''
Returns a launcher that when run on the local machine
will transfer input dicts appropriately and execute the 
current file with an input parameter specifying arguments.

'''

    scriptfile = os.path.abspath(inspect.stack()[0][1])
    scriptroot = 'compbio'
    func = 'remote_batch_1d'
    launcher = bsub.local_launcher(scriptfile,
                                   scriptroot,
                                   func = func,
                                   input_dicts = input_dicts,
                                   run_id = run_id,
                                   host = host)
    return launcher
                       

#ALT Remote launchpoint for bsub.
def remote_batch_1d(run_id):
    '''
Remote job mothership for a single run that will generate
a bunch of bsub commands to run a series of RNAz checks
along the gene of interest.

inputs:
  run_id

output:
  the datapath (same for local and remote) of data output each thread
'''

    inp_dicts = butils.load_data(run_id, 'input')    
    eyeball = bsub.eyeball(run_id,
                           os.path.abspath(inspect.stack()[0][1]), inp_dicts,
                           func = 'remote_run_1d',
                           name = run_id+'_eye', 
                           mem = 2)
                          
    eyeball.launch()
    eyeball.await()
    eyeball.package()
    eyeball.complete()

#Remote launchpoint for bsub.
def remote_run_1d(run_id):
    '''
    

'''
    baserng = p['baserng']
    win_len = p['win_len']
    win_ofs = p['win_ofs']
    
    ali_nums = day0.fetch_num_ali()
    l = shape(ali_nums[0])
    

    


        
def usage():
  print '''
usage: 
day.py function run_id

Call function with run_id.
'''
  exit(1)

if __name__ == '__main__':
    run_id = sys.argv[2]
    run_func = globals()[sys.argv[1]]
    output_dict = run_func(run_id)
    if output_dict == None:
        output_dict = {'blank':'Nothing output in call to {0}'.\
                           format(sys.argv[1])}
    butils.save_data( output_dict, run_id, 'output')
