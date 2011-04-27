================
tophatBedToJuncs
================
Takes a Tophat junctions.bed file and converts it to the raw .juncs 
format, adding a junction-location unique junction ID. This makes for
easy comparison of junctions across different junction files using
standard command line utilities.

======================
tophatBedToSeq
======================
Takes a Tophat junctions.bed file and converts it to a FASTA file of
sequences specifying the junctions. Takes a parameter to determine
the size of the sequences to return. It is required that fastaFromBed
from the BEDTools_ package is installed and in your path.

=================
filterCuffCompare
=================
Takes a tracking file from Cuffcompare and a combined.gtf file from
Cuffcompare and retains transcripts with a specified number of samples
supporting evidence of the transcript. This also retains all transcripts
that match the reference annotation. Useful for wiping out a lot of the
junky partially-assembled transcripts that Cufflinks spits out.

.. _BEDTools: http://code.google.com/p/bedtools/
