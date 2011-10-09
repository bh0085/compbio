import cb.utils.memo as mem
import cb.p.network.io as nio
import cb.p.nfilt.comparator as comp


def align_graphm():
    nets = comp.get_graphs(restriction = 'bdtnp')    
    return nets
