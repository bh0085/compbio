'''
Load up the part of the 32 mammals genome corresponding to PvT1.

Check coverage.

Save a CLW file.

Run RNAz on windows.

Then run RNAz for long range interactions.

'''
import cb.config as cfg
import cb.utils.plots as myplots
import Bio.AlignIO as aio

def run0():
    bases = (128693265,129266680)
    fa_file = cfg.dataPath('pvt1/pvt1.fa')

    
    ali = aio.parse(open(f), 'fasta')
    a0 = ali.next()

    f = myplots.fignum(3,(8,8))
    ax = f.add_subplot(111)
    vec = zeros(len(ax[0]))
    raise Exception()

    f.savefig(myplots.figpath('run0_alignment_hits'))

    
