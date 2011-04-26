'''remote_utils.py

A bunch of scripts for doing things on the remote machine including:

ssh_exec()
  Take a shell command, quote it appropriately and execute it via 
  ssh.

ssh_cmd()
  Taks a shell command, quote it and prepend ssh for remote exec
  via subprocess.Popen.

remote_datapath():
  Uses the config.dataPath() function on remote and/or local machines
  to determine correct absolute paths for a datafile. Output can be
  used to pinpoint files so that they can be sshd as in 'scp_data()'.

  It is called with a timeout and so if the datapath query takes too
  long to execute (presumably because ssh authentication has failed)
  raises an exception.

scp_data():
  scps a file from src to dest with paths given relative to the 
  /data directory on the src/host machines specified.
  
  Intended usage is to corral files from remotely executed scripts
  when script output dictionaries contain an 'outfile' keyword.

  See bsub_utils.py class local_launcher::fetch_start()

comm_timeout():
  Execute a subprocess with a timeout timer.
  Useful for ssh when authentication may fail and script will
  hang.
'''

import pipes , subprocess as spc
import signal
from signal import alarm, SIGALRM, SIGKILL
import compbio.config
from os import kill


def ssh_exec( cmd , host = None):
    sshcmd = ssh_cmd(cmd, host = host)
    prc = spc.Popen(sshcmd, shell = True, stdout = spc.PIPE)
    return prc.communicate()[0]
    


def ssh_cmd(cmd, host = None):
    if host == None:
        return cmd
    else:
        return  '''ssh {0} {1}'''.\
            format(host, pipes.quote(cmd))



def remote_datapath( path, host, volume = '' , timeout = 100):
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
    return prc




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



class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm

def get_process_children(pid):
    p = spc.Popen('ps --no-headers -o pid --ppid %d' % pid, shell = True,
              stdout = spc.PIPE, stderr = spc.PIPE)
    stdout, stderr = p.communicate()
    return [int(p) for p in stdout.split()]
