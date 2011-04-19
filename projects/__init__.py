#!/usr/bin/env python
'''
Projects module. 

Here lives the compbio projects thatI  have under active construction that are not ready to graduate. To freeliving modules.

Nearly ready to graduate are:

  cbdb: A tool for storing biological data in sqlite database mapped to sqlalchem classes. Various utilities and helper functions are available custom tailored to particular large databases e.g: ncbi taxonomy, genbank sequence files.

Nascent but having substantial work in them are:

  seqtree: A tree browser that reads tree in from a 16s tree of microbial life and has utilities to grab sequences out of alignments at taxonomical nodes as specified in the join table from ncbi (ncbi_gbjoin in cbdb). Seqtree has utilities to grab sequences belonging to a given ncbi clade, to create per clade alignments, and eventually, to gauge varying evolution rates for example for a given ghene family compaored to the 16s mutation rate as they vary across edges int eh tree of life.

  GC: A stub that so far mostly demonstrates a use case for the genome browser I have been debeloping in draw_hilbert. Currently, plots GC skew for a series of genomes in a hilbert curve/cloud format from genbank files. There is much room to build. At some point, I hope to incorporate mutation repair and thoughts about codon bias etc. I am also interested in relating this to recent development in rotating proteins. Is it possible that a prescence and abscence of DNA bindign rotatry proteins correlates to mutation repair efficiency? Can GC skew be regarded as a proxy for DNA repair efficiency?

  predict: My gene predictor. Contains some generic classes for machine lerning. I am interested now in putting in a genetic algorithm for novel edge detection in a gene regulatory network such as patrick's. Have some notes tied to a discussion with dave living in keepnote.

'''
