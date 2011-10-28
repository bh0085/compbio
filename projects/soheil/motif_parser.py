#!/usr/bin/env python
import itertools as it

def parse(file='/data/soheil/regions.txt'):
    fopen = open(file)
    motifs = []
    while(1):
        l1 = fopen.readline()
        #l2 = fopen.readline()

        if l1.strip() == '': break;
        fields = [e.strip() for e in l1.strip().split('\t')]
        seq = fields[7]
        score = int(fields[6]) if fields[6] != '-' else -1
        
        motifs.append({'score':score,
                       'name': fields[5],
                       'seq':seq,
                       'line': l1})

    return motifs

def best_30(motifs):
    factor_lines = dict([(k,list(g)) 
                         for k, g in it.groupby(sorted(motifs,key = lambda x:x['name']),
                                                key = lambda x: x['name'])
                         ])

    outputs = {}
    for fname, lines in factor_lines.iteritems():
        flt_motifs = [m for m in lines if m['score'] != -1 \
                          and len(m['seq']) < 2000]
        srt_motifs = sorted(flt_motifs, key = lambda x: x['score'])
        outputs[fname] = srt_motifs[:30]

    return outputs

def usage():
    print '''usage:
./motif_parser file1 [file2 file3...]

Reads in the file arguments and outputs best motifs for every
factor in the dataset.

one way to use this would be:

ls *txt | xargs ./motif_parser

to rapidly process all text files in a directory.
'''
if __name__ == '__main__':
    import sys
    import os
    args = sys.argv
    files = args[1:]
    if len(files) == 0: exit(usage())
    for f in files:
        motifs = parse(f)
        outputs = best_30(motifs)

        for fname, matches in outputs.iteritems():
            lines = [b['seq'] for b in matches]
            new_fname = '.'.join(list(os.path.splitext(f)[:-1]) + \
                                     ['{0}.top30.txt'.format(fname)])
            fopen = open(new_fname, 'w')
            fopen.write('\n'.join([line for line in lines]))
            fopen.close()

        print 'done with {0}'.format(f)
        
