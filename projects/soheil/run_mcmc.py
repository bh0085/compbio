#!/usr/bin/env python
import sys
import compbio.utils.bsub_utils as bsu

def launch_many(run_id):
    '''
Generate script paramaters and launch a bunch of bsub jobs.

Designed to be run on the cluster via an interactive shell.
Note: If this is not run on cluster, since it does not look
up a remote url for files, it won't be able to find expression
data.

'''
    print 'Launching all jobs!'

    #MAKE INPUTS 
    expr_filename = 'soheil/expression_c4d_n4_tt_{0}.mat'.format(ttnum)
    url = cfg.dataURL(expr_filename)
    remote_exprname = cfg.dataPath(url)

    inp_dicts = []
    for iters in [2,4,8]:
        inp_dict = dict(out_iter_num = iters,
                        in_iter_num = 5,
                        k = 6,
                        beta = 4,
                        f_mix = 1,
                        f_sim = .9,
                        f_en_in = .95,
                        f_en_out = .95,
                        th_cor = .5,
                        trunc_value = 3,
                        degree_bound = 3,
                        filename = remote_exprname)
        inp_dicts.append(inp_dict)
    

    #MAKE EYEBALL
    eyeball = bsub.eyeball(run_id, 
                           os.path.abspath(inspect.stack()[0][1]),
                           inp_dicts,
                           func = 'run_single',
                           name = 'mcmc_',
                           mem = 3)

    #LAUNCH EYEBALL JOBS
    eyeball.launch()

    
    #RETURN A LIST OF LAUNCHED JOBS
    return dict(cmds=eyeball.cmds)



def run_single(run_id):
    '''
Given an input dictionary containing a single paramater set
run mcmc in matlab using bs_macros.run_matlab.
'''
    input_dict = butils.load_data(run_id,'input')
    return bsm.runmat('run_mcmc', input_dict, run_id)




if __name__ == '__main__':
    run_id = sys.argv[2]
    run_func = globals()[sys.argv[1]]
    output_dict = run_func(run_id)
    if output_dict == None:
        output_dict = {'blank':'Nothing output in call to {0}'.\
                           format(sys.argv[1])}
    bsu.save_data( output_dict, run_id, 'output')
