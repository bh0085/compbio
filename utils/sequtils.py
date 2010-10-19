
import re
n_lookup={'A':'A',
        'T':'T',
        'G':'G',
        'C':'C',
        'U':'U',
        'R':'[AG]',#Purine
        'Y':'[CT]',#Pyramidine
        'N':'[AGCT]',#Nucleotide
        'W':'[AT]',#Weak
        'S':'[GC]',#Strong
        'M':'[AC]',#Amino
        'K':'[GT]',#Keto
        'B':'[GCT]',
        'H':'[ACT]',
        'D':'[AGT]',
        'V':'[AGC]'
        }

lookup_degeneracy ={
    'A':1,
    'T':1,
    'G':1,
    'C':1,
    'U':1,
    'R':2,#Purine
    'Y':2,#Pyramidine
    'N':4,#Nucleotide
    'W':2,#Weak
    'S':2,#Strong
    'M':2,#Amino
    'K':2,#Keto
    'B':3,
    'H':3,
    'D':3,
    'V':3
    }


def sequence_re(seq):
    #assume that motifs represents a list of motif possibilities
    #with wild cards.
    r = ''
    upper = seq.upper()
    for s in upper:
        r += n_lookup[s]
    import re
    return re.compile(r,re.I)

    
def kmer_in_list(kmer,motifs):
    k = len(kmer)

    max_score = 0.0
    max_index = -1
    for i in range(len(motifs)):
        m = motifs[i]
        
        match_score = 0.0
        match_idx = -1 
        if k <= len(m):
            for j in range( len(m) - k + 1):
                sub = m[j:j+k].upper()
                r = sequence_re(m[j:j+k])
                if re.search(r,kmer):
                    s = 0.0
                    for letter in sub: s += 1.0/lookup_degeneracy[letter]
                    if s > match_score:
                        match_idx = j
                        match_score = s
                    
        else:
            r = sequence_re(m)
            if re.search(r, kmer):
                for letter in m: s += 1.0/lookup_degeneracy[letter]
                match_idx = 0 
        if match_idx != -1 and match_score > max_score:
            max_score = match_score
            max_index = i

    return (max_index, max_score)
