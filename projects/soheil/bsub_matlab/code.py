#!/usr/bin/env python
'''
Code to run a bunch of matlab scripts
on the cluster using parameters generated
(by hardcoded routines) in python.
'''
#IMPORTS
import scipy.io as sio
import sys
import subprocess 

#GLOBALS
#ALSO MODIFY THIS PART-------------------------------
run_type = 'TEST'

def get_stuff_for_run():
    '''
usage:
   
infiles,outfiles,params,run_ids = get_stuff_for_run()

'''
    #Put code in here to generate the 
    #List of input and output files in your
    #directory system as well as the params 
    #that will be saved in the input files
    #and eventually fed to your matlab script.

    if run_type == 'TEST':
        nruns = 2
        infiles = ['inp_test_{0}'.format(i) for i in range(nruns)]
        outfiles = ['out_test_{0}'.format(i) for i in range(nruns)]
        params = [{'number': i,
                   'motif_file': 'some_motif_file.txt'} 
                  for i in range(nruns)]
        run_ids = ['TEST{0}'.format(i) for i in range(nruns)]



    #MODIFY THIS PART ONLY!---------------------------
    elif run_type == 'YOURRUNTYPEHERE':
        '''just make sure that you set values for the four
        crucial parameters:
        infiles, outfiles, params, run_ids'''
        

        pass
    #END MODIFIED PART HERE!---------------------------


    else:
        print 'RUN TYPE NOT UNDERSTOOD'
    return infiles, outfiles, params, run_ids




def run(script, input_file, output_file, params,  run_id):
    '''
usage:
run(input_file, output_file, params)

inputs:
  filenames:   [input_file], [output_file]
  script name: [whaveter]
'''

    #SAVE PARAMS TO INP FILE
    sio.savemat(input_file,params)
    
    #THIS IS THE COMMAND THAT MATLAB WILL RUN
    mat_cmd = '''\\"{2}('{0}', '{1}' ); exit\\"'''.\
        format(input_file,output_file,  script)

    #THIS IS THE COMMAND THAT WILL RUN MATLAB
    cstr = '''echo {0} | matlab -nojvm -nodisplay -nosplash '''.\
        format(mat_cmd)

    #THIS IS THE COMMAND THAT THE CLUSTER WILL RUN 
    #TO RUN THE SHELL COMMAND THAT RUNS MATLAB
    bsc = 'bsub -o {0} -q compbio-week -P {1} "{2}"'.format(\
        'logfile_{0}'.format(run_id), script, cstr)

    #THIS THE COMMAND WITHIN PYTHON
    #THAT LAUNCHES THE CLUSTER CMD
    #....
    subprocess.Popen(bsc, shell = True)




if __name__ == '__main__':
    script = sys.argv[1]
    infiles, outfiles, all_params, run_ids = get_stuff_for_run()
    for args in zip( infiles, outfiles, all_params,run_ids):
        run(script, *args)
