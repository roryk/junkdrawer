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
partially-assembled transcripts that Cufflinks spits out. Takes
-s, the number of samples a transcript must exist in, -f a gtf file and
-t a tracking file.

The second mode filters transcripts based on their proportion of
exons in the longest transcript for that locus. Takes -f a gtf file and
-e the threshold proportion. For example a locus with 4 exons will
pass if 3 are in the transcript and fail if only 1 exon is in the
transcript if -t is set to .30.

editAttributes
==============
Deletes or swaps attributes in the attribute line of a GTF file. -d deletes the
attributes. -r and -s replaces the attributes in -r with the ones in -s.

.. _BEDTools: http://code.google.com/p/bedtools/
