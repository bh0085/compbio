#!/usr/bin/env python

'''
bsjobs.py

Scripts designed to be called from the command line.
With jobids as input.
'''

from compbio.utils.bsub_utils import *

def bjobs(job_dict):
    '''
Return the run statuses of programs having given jobids. Uses bjobs.
'''

    if len(job_dict) == 0:
        return {}
    jobids  =job_dict.values()
    jobnames=job_dict.keys()
    job_idnames=dict([(v,k) for k,v in job_dict.iteritems()]) 
    job_stats = {}




    #Get active jobs.
    if jobids.count(-1) != len(jobids):
        jobs = subprocess.Popen('bjobs '+ ' '.join(['{0}'.format(j) for j in jobids
                                                    if not j == -1]), 
                                shell = True, 
                                stdout = subprocess.PIPE).\
                                communicate()[0]
        lines = jobs.split('\n')
        cols, lines = lines[0],lines[1:]
        col_starts = {}

    #Parse job list into columns
        terms =  re.compile('\s+').split(cols)
        starts = [cols.index(t) for t in terms]
        ends =   [cols.index(t) for t in roll(terms, -1)]
        ends[-1] = len(cols)
        col_ranges = dict([(t, (starts[i], ends[i])) 
                       for i, t in enumerate(terms)])
    
        job_dicts = {}
        for l in lines:
            if l.strip() == '': continue
            d0 = dict([(k, l[col_ranges[k][0]:col_ranges[k][1]]) 
                       for k in terms])
            run_id = int(d0['JOBID'])
            job_stats[job_idnames[run_id]] = d0

      
    for k in jobnames:
        if not k in job_stats.keys():
            job_stats[k] = {'STAT':'UNSUBMITTED',
                            'JOBID':job_dict[k]}


    return(job_stats)


if __name__ == '__main__':
    assert len(sys.argv) > 2
    if sys.argv[1] in [ 'bjobs']:
        #Dump data to json and write to stdout
        ids = sys.argv[2:]
        sys.stdout.write(sjson.dumps(globals()[sys.argv[1]](ids)))
        exit(0)
    else:
        raise Exception()
