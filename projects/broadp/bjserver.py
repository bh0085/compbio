import pp
import subprocess
import itertools as it
import time

nodes =[]
class BJServer():
    def __init__(self, name = 'BJDefault', \
                     worker_timeout = 10, \
                     secret_key = 'secret_key',\
                     n_workers = 2):
        self.job_server = None
        self.jobs = []
        self.nodes = []
        self.name = name
        self.secret_key = secret_key
        self.worker_timeout = worker_timeout
        self.n_workers = n_workers
        self.initialized = False

    def begin(self):
        self.launch_workers()
        nodes, ids = self.await_workers()
        self.make_manager(nodes, ids)
        
    def end(self):
        self.kill_workers()
        self.job_server.destroy()

    def launch_workers(self):
        number = self.n_workers
        worker_task = '''bsub -P %s -o bsub_out/outfile$i -q compbio-week "ppserver.py -t %s -a"''' % (self.name, self.worker_timeout)
        command = "ssh tin 'for i in {1..%s}; do %s ; done'" % \
            (number, worker_task)
        print command
        prc = subprocess.Popen(command, shell = True, \
                                   stdout = subprocess.PIPE)
        return self.await_workers()

    def check_workers(self):
        command = "ssh tin 'bjobs -P %s'" %self.name
        print command
        prc = subprocess.Popen(command,shell = True, 
                               stdout = subprocess.PIPE)
        output = prc.communicate()[0]
        lines = output.split('\n')
        if lines == ['']: return [], []
        name_rng = (lines[0].index('EXEC_HOST'),
                    lines[0].index('JOB_NAME'))
        id_rng = (lines[0].index('JOBID'),
                  lines[0].index('USER'))

        nodes_running = map(lambda x: x[name_rng[0]:name_rng[1]],
                            (it.ifilter(lambda x : 'RUN' in x,lines[1:]))
                            )
        ids_running = map(lambda x: x[id_rng[0]:id_rng[1]],
                            (it.ifilter(lambda x : 'RUN' in x,lines[1:]))
                            )
        print lines
        return nodes_running, ids_running
    def await_workers(self):
        start_time = time.time()
        while 1:
            if time.time() > start_time+100:
                raise Exception('timed out... workers failed to come online')

            nodes, ids = self.check_workers()
            if len(nodes) < self.n_workers:
                time.sleep(.5)
                continue
            break

        print 'success, all workers online at nodes:'
        for n in nodes: print n
        return nodes, ids
        
    def make_manager(self, nodes, ids):
        assert not self.initialized
        servers = nodes
        servers = tuple(map(lambda x: str.strip(x), servers))
        job_server = pp.Server(ppservers = servers)

        these_nodes = [{'name':nodes[i],
                       'id':ids[i],
                       'manager':job_server}\
                          for i in range(len(nodes))]
    
        globals()['nodes'].extend(these_nodes)
        self.nodes = these_nodes


        self.job_server = job_server
        self.initialized = True

    def kill_workers(self):
        assert self.initialized
        kill_ids = []
        for j in self.nodes[::-1]:
            print j['manager']
            if j['manager'] == self.job_server:
                kill_ids.append(j['id'])
                self.nodes.remove(j)
        kill_cmd = "ssh tin 'for id in %s ; do  echo ${id}; bkill ${id}; done'"\
             % ' '.join(kill_ids)

        print kill_cmd
        prc = subprocess.Popen(kill_cmd, shell = True, stdout = subprocess.PIPE)

        comm = prc.communicate()
        print comm
        self.initialized = False
    
    def integrate_cos(k):
        nx = 1000
        f1 = numpy.cos(arange(nx))
        f2 = numpy.cos(arange(nx)*k)
        val = numpy.sum(f1 * f2)
        return val

    def submit(self,func, name = None, **kwargs):
        if not name:
            name = func.__str__() + kwargs.get('args', ()).__str__()
        self.jobs['name'] = self.job_server.submit(func, **kwargs)
