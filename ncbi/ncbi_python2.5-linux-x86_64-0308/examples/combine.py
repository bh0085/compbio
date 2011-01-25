#!/usr/bin/env python

# This script takes a set of text asn.1 Seq-entry files,
# checks whether they're nucleic acids, combines those that are
# into a Bioseq-set, wraps this in a Seq-entry, adds the date
# and a comment, and writes it out as asn.1 text.
#
# Try, for example:
#     combine.py data/*.prt > my_entry.prt

import ncbi
import sys, os

seq_set = ncbi.CBioseq_set()

files_used = []
for input_file in sys.argv[1:]:
    entry = ncbi.ReadAsnFile(input_file)
    id = entry.GetSeq().GetFirstId().AsFastaString()
    if entry.GetSeq().IsNa():
        seq_set.SetSeq_set().append(entry)
        sys.stderr.write('added %s from file %s\n' % (id, input_file))
        files_used.append(input_file)
    else:
        sys.stderr.write(
            'WARNING: skipping %s from file %s; '
            'not a nucleic acid sequence\n'
            % (id, input_file))

seq_set.SetDate(ncbi.CDate(ncbi.CurrentTime()))

comment = ncbi.CSeqdesc()
comment.SetComment('Produced by user %s using combine.py on files %s'
                   % (os.environ['LOGNAME'], files_used))
seq_set.SetDescr().Set().append(comment)

big_entry = ncbi.CSeq_entry()
big_entry.SetSet(seq_set)

ncbi.cout << ncbi.MSerial_AsnText << big_entry;
