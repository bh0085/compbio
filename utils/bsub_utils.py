#/usr/bin/env python

import sys  

def bjobs(jobids):
    '''
    Return the run statuses of programs launched under the control of
this eye. Uses bjobs.
'''
    #Get active jobs.
    jobs = subprocess.Popen('bjobs '+ ' '.join(jobids), 
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
      run_id = d0['JOBID']
      job_dicts[run_id] = d0

if __name__ == '__main__':
    print 'hello'
    exit(0)
    assert len(sys.argv) > 2
    if sys.argv[1] == 'bjobs':
        ids = sys.argv[2:]
        print ids
        exit(0)
    else:
        raise Exception()
