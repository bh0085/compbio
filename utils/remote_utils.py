import pipes , subprocess as spc
import signal
from signal import alarm, SIGALRM, SIGKILL

from os import kill


def ssh_exec( cmd , host = 'tin'):
    sshcmd = ssh_cmd(cmd, host = host)
    prc = spc.Popen(sshcmd, shell = True, stdout = spc.PIPE)
    return prc.communicate()[0]
    
def ssh_cmd(cmd, host = 'tin'):
    return  '''ssh {0} {1}'''.\
        format(host, pipes.quote(cmd))

def scp_data(src_url, dest_url, timeout = 5):
    prc = spc.Popen('scp -v {0} {1}'.\
                        format(cfg.dataPath(src_url),
                               cfg.dataPath(dest_url)),
                    shell = True)
    result = prc.communicate()

def remote_datapath( path, host ):
    timeout = 5

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
