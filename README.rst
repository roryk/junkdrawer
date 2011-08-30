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
filter on an attribute and reorder the attributes field.

ensembl_cleaner
===============
Fixes the GTF file from Ensembl for use with Tophat and Cufflinks by putting
in a consistent set of chromosome names.

transcriptLength
================
Calculates the length of all transcripts in a GTF file by adding up all
of the 'exons' features for each transcript. Also can filter a GTF file
by minimum and maximum transcript lengths. This is useful for removing
small transcripts that may be lost during RNA purification such as 
microRNA and snoRNA and also to remove large transcripts that may be
degraded.

transcript_per_gene_counter
===========================
counts up the number of transcripts per gene in a GTF file. the GTF file
must have the format gene_id foo; transcript_id bar, otherwise it won't
work.

gtf2bed
=======
converts a gtf file to a bed file of unique intervals

gtf2fasta
=========
converts a gtf file to a file of fasta sequences of the transcripts.
it is required that the attribute list contains a transcript_id.

combineJuncs
============
combines identical junctions, adding their score together. Workflow to
use this is to convert the Tophat BED files to juncs with tophatBedToJuncs
and then combine the .juncs files with combineJuncs.

juncs2seq
=========
converts a juncs file to a file of sequences for the junctions. Usage
is juncs2seq pathtofastafile basestokeep< juncsfile. Bases to keep
is the bases to keep on each side of the splice junction. Requires
samtools to be installed and in your path. Skips junctions that go
through repeat regions (regions with lowercase letters in the genome 
file).

gtf2juncs
=========
converts a GTF file to a juncs file, outputting all possible
single exon skipping events in and all splice junctions in 
juncs format.


.. _BEDTools: http://code.google.com/p/bedtools/
