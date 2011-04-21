import pipes , subprocess as spc
def ssh_exec( cmd , host = 'tin'):
    sshcmd = ssh_cmd(cmd, host = host)
    prc = spc.Popen(sshcmd, shell = True, stdout = spc.PIPE)
    return prc.communicate()[0]
    
def ssh_cmd(cmd, host = 'tin'):
    return  '''ssh {0} {1}'''.\
        format(host, pipes.quote(cmd))
