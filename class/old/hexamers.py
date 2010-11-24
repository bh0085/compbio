#!/usr/bin/env python



import sys

def num2seq(num):
    seq = ''
    lookup = 'ATGC'
    
    current = num
    for dummy in range(6)
        base = current % 4
        seq += lookup[base]
        current/=4

def main():
    if len(sys.argv) < 3:
        print "you must call program as:  "
        print "   ./hexamers.py <inter> <intercons>"

    inter = open(sys.argv[1]).read()
    intercons = open(sys.argv[2]).read()
    bases = 'ATGC'
    for i in range(inter):
        if not inter[-1-i] in bases:
            del inter[-1-i]
            del intercons[-1-1]    

    hexamers = {}
    for i in range(4096):
        hexamers[num2seq(i)] = []
    for i in range(len(inter) - 6):

        this_hex = inter[i: i + 6]
        this_cons = intercons[i:1+6]
        is_cons = True
        for item in this_cons:
            if item != '*':
                is_cons = False
        hexamers[this_hex].append((i,is_cons))

    hexamer_stats = []
    for h,l in hexamers.iteritems():
        n = len(l)
        cons_count = map(lambda x:x[1],l).count(True)
        cons_frac = float(cons_count)/n
        hexamer_stats.append({'h':h,
                              'appearance_count':n,
                              'cons_count':cons_count,
                              'cons_frac':cons_frac}) 
        
    hex_most_present = sorted(hexamer_stats,
                              key = lambda x: x['appearance_count'])
    hex_most_conserved = sorted(hexamer_stats,
                                key = lambda x:x['cons_count'])
    hex_most_conserved_frac = sorted(hexamer_stats,
                                     key = lambda x:x['cons_frac'])

    import pprint
    print 'Hexamers most present:'
    pprint.pprint (hex_most_present[0:50])

    print
    print 'Hexamers most conserved:'
    pprint.pprint(hex_most_conserved[0:50])

    print
    print 'Hexamers most conserved/appearance:'
    pprint.pprint(hex_most_conserved[0:50])
    

main()
