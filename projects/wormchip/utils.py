import cb.config as cfg
import cb.utils.memo as mem

def get_map_rows(**kwargs):
    def set_map_rows(**kwargs):
        mapfile = cfg.dataPath('wormbase/loci_all.txt')
        fopen = open(mapfile)
        lines = fopen.readlines()
        cols = [e.strip() for e in lines[0].strip().split(',')]
        rows =[ dict( zip(cols,[e.strip() for e in l.strip().split(',')])) 
                for l in lines[1:-1]]
    
        return rows
    return mem.getOrSet(set_map_rows, **mem.rc(kwargs))

def symbol_ids():
    rows = get_map_rows()
    return dict([(r['Gene name'],r['WormBase gene ID']) for r in rows])

