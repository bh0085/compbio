import os
#Wraps unafold.
def fold(seq):

    #Calls UNAFold's hybrid folding routine and returns dG
    #in addition to the sequence matches.

    
    if not os.path.isdir('temp'):
        os.path.mkdir('temp')

    os.chdir('temp')
    f = open('unatemp.seq','w')
    f.write(seq)
    f.close()
    import subprocess
    subprocess.call(["hybrid-ss-min unatemp.seq"],shell=True)
    fout = open('unatemp.ct')

    os.chdir('..')

    header = fout.readline()
    rows = fout.read()
    import re
    parsed = re.finditer(re.compile('\w+\t\w+\t\w+\t\w+\t(\w+).*',re.M),rows)
    found = list(parsed)
    pairs = []
    
    dg = re.search(re.compile('dG = ([^\t]*)'),header).group(1)
    
    for item in found:
        pairs.append(int(item.group(1)))
        
    return pairs, dg
