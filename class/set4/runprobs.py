import tree
import wfsim
import snp

def run(name = 'tree'):
    if name == 'tree':
        tree.run()
    elif name == 'wfsim':
        wfsim.run()
    elif name == 'snp':
        snp.run():
    else :
        raise Exception("unhanded program type")


