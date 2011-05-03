tophatBedToJuncs
================
Takes a Tophat junctions.bed file and converts it to the raw .juncs 
format, adding a junction-location unique junction ID. This makes for
easy comparison of junctions across different junction files using
standard command line utilities.

tophatBedToSeq
======================
Takes a Tophat junctions.bed file and converts it to a FASTA file of
sequences specifying the junctions. Takes a parameter to determine
the size of the sequences to return. It is required that fastaFromBed
from the BEDTools_ package is installed and in your path.

filterCuffCompare
=================

Two modes, the first takes a tracking file from Cuffcompare and a
combined.gtf file from Cuffcompare and retains transcripts with a
specified number of samples supporting evidence of the
transcript. This also retains all transcripts that match the reference
annotation. Useful for wiping out a lot of the junky
partially-assembled transcripts that Cufflinks spits out. 

The second mode filters transcripts based on their proportion of
exons in the longest transcript for that locus. Takes a gtf file and
the threshold proportion. For example a locus with 4 exons will
pass if 3 are in the transcript and fail if only 1 exon is in the
transcript if the threshold is set to .30.

editAttributes
==============
Deletes, swaps, adds GTFs to the attribute line of a GTF file. Can also
filter on an attribute.

ensembl_cleaner
===============
Fixes the GTF file from Ensembl for use with Tophat and Cufflinks by putting
in a consistent set of chromosome names.

.. _BEDTools: http://code.google.com/p/bedtools/
