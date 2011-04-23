import pipes , subprocess as spc
import signal
from signal import alarm, SIGALRM, SIGKILL
import compbio.config
from os import kill


def ssh_exec( cmd , host = 'tin'):
    sshcmd = ssh_cmd(cmd, host = host)
    prc = spc.Popen(sshcmd, shell = True, stdout = spc.PIPE)
    return prc.communicate()[0]
    
def ssh_cmd(cmd, host = 'tin'):
    return  '''ssh {0} {1}'''.\
        format(host, pipes.quote(cmd))

def remote_datapath( path, host, volume = '' ):
    timeout = 5
    prc = spc.Popen('ssh {0} "pydatapath.py {1} {2}"'.\
                        format(host, path, volume), 
                    shell = True, stdout = spc.PIPE)    
    out = comm_timeout(prc, timeout)
    if out[0] == -9:
        raise Exception('Process timeout')
    else: 
        return out[1]

def scp_data(src_path, dest_path,
             src_host = None, dest_host = None):
    cfg = compbio.config
    if src_host:
        src_url = src_host + ':' +cfg.dataPath(\
            cfg.dataURL(src_path, host = src_host))
    else:
        src_url = cfg.dataPath(src_path)
    if dest_host:
        dest_url = dest_host + ':' + cfg.dataPath(\
            cfg.dataURL(dest_path, host =dest_host))
    else:
        dest_url = cfg.dataPath(dest_path)

    prc = spc.Popen('scp {0} {1} '.\
                        format(src_url, dest_url), 
                    shell = True, stdout = spc.PIPE)    
    return prc.communicate()

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm

def comm_timeout(proc, time):
    signal.signal(SIGALRM, alarm_handler)
    signal.alarm(time)  # 5 minutes
    try:
        stdoutdata, stderrdata = proc.communicate()
        signal.alarm(0)  # reset the alarm
        return proc.returncode, stdoutdata, stderrdata
    except Alarm:
        pids = [proc.pid]
        kill_tree = 1
        if kill_tree:
            pids.extend(get_process_children(proc.pid))
        for pid in pids:
            kill(pid, SIGKILL)
        return -9, '', ''

def get_process_children(pid):
    p = spc.Popen('ps --no-headers -o pid --ppid %d' % pid, shell = True,
              stdout = spc.PIPE, stderr = spc.PIPE)
    stdout, stderr = p.communicate()
    return [int(p) for p in stdout.split()]
