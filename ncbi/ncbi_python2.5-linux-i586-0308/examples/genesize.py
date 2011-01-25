# Example script for the NCBI C++ Toolkit Python interface.
# Run by typing execfile('genesize.py') at the python prompt,
# or python genesize.py at the unix prompt.
#
# For each human chromosome, iterate over the annotated
# genes, recording their total lengths.  Report the
# mean and maximum lengths for each chromosome.
# Show the largest gene for each chromosome in
# a web browser.

import ncbi

scope = ncbi.CSimpleOM.NewScope()

lengths = []; mean_lengths = []; max_lengths = []; largest_genes = []
counts = []; names = []
print '%s\t%s\t%s\t%s\t%s' % ('\nChrom.', 'mean', 'max', 'max id', 'genes')
for i in range(1, 25):
    # Use the fact that the human chromosomes have accessions of the
    # form 'nc_0000dd', where 'dd' ranges from 01 to 24.
    acc = 'nc_0000%.2d' % i
    name = str(i)
    if i == 23: name = 'X'
    if i == 24: name = 'Y'
    id = ncbi.CSeq_id(acc)
    hand = scope.GetBioseqHandle(id)

    sel = ncbi.SAnnotSelector(ncbi.CSeqFeatData.e_Gene)
    sel.SetResolveAll().SetResolveDepth(1)

    tmp_lengths = []; max_length = 0; gene_count = 0
    gene = ncbi.CFeat_CI(hand, sel)
    while gene:
        loc = gene.GetLocation()
        length = loc.GetInt().GetTo() - loc.GetInt().GetFrom() + 1
        tmp_lengths.append(length)
        if length > max_length:
            max_length = length
            largest_gene = gene.GetDbxref()[0].GetTag().GetId()
        gene_count += 1
        gene.Incr()
    lengths.append(tmp_lengths)
    counts.append(gene_count)
    mean_length = float(reduce(lambda x, y: x + y, tmp_lengths)) \
                  / len(tmp_lengths)
    mean_lengths.append(mean_length)
    max_lengths.append(max_length)
    largest_genes.append(largest_gene)
    names.append(name)
    print '%s\t%.0f\t%d\t%d\t%d' \
          % (name, mean_length, max_length, largest_gene, gene_count)

# print 'em sorted by mean gene length
print '%s\t%s\t%s\t%s\t%s' % ('\nChrom.', 'mean', 'max', 'max id', 'genes')
idx = range(0, len(mean_lengths))
idx.sort(lambda l, r: cmp(mean_lengths[r], mean_lengths[l]))
for i in idx:
    print '%s\t%.0f\t%d\t%d\t%d' % \
          (names[i], mean_lengths[i], max_lengths[i], \
           largest_genes[i], len(lengths[i]))

# launch a web page showing the longest gene for each chromosome
ncbi.EntrezWeb(largest_genes, 'gene')
