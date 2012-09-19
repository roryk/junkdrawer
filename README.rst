Example workflows
=================
The experiment: you have a single lane of data, and have gone through
the Tophat process.

To find all splice junctions that are novel to an annotation, assuming
you have a GTF file for the annotation:

1. tbed2juncs < junctions.bed > junctions.juncs
2. gtf2juncs < your_annotation.gtf > annotated.juncs
3. cat junctions.juncs | cut -f4 | sort | comm -23 - annotated.juncs >
novel_junctions.txt

To find all splice junctions that map to an EST database that you
made with bowtie, allowing for repeat regions.
1. juncs2seq rn4.fa 25 0 < junctions.juncs > junctions.fa
2. bowtie -v 2 rat_ests junctions.fa | cut -f1 > hit_ests.txt

To enumerate all splice junctions from a GTF file including all
possible exon skipping events with the canonical splicing signals
GC-AG, AT-AC, GC-AG :
1. gtf2juncs --skip < annotation.gtf | uniq | sort | uniq > juncs.juncs
2. juncs2seq genome.fa -2 --i --rt < juncs.juncs > juncs.fa
3. egrep -i -B 1 -e "GTAG|ATAC|GCAG" juncs.fa | grep -ve '--' > canonical.fa

To enumerate all possible inclusion/exclusion events for an GTF file and
find counts for each event from a juncs file:
1. gtf2juncs --ie < annotation.gtf > annotation.ie
2. ie_calculator.py annotation.ie juncs.juncs > ie.counts

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

tbed2juncs
==========
convert tophat bed file to a juncs file

grepGTF
=======
grep a GTF file for a field or attribute with IDs from another file.
For example to find all entries in the file ensembl.gtf that have the gene_ids
listed in the file tofilter.tsv:

python grepGTF ensembl.gtf gene_id tofilter.tsv

To find the ones that don't match add -v to th end.

shuffle_fasta
=============
Shuffles the sequences in a FASTA file, good for quick and dirty false-positive tests.
Takes a FASTA file as stdin and returns shuffled sequences, each sequence maintaining
the same overall nucleotide composition.

detect_fastq_format
===================
Outputs a list of possible formats a FASTQ file could be after checking at most
a million reads.

usage: detect_fastq_format.py fastq_file

.. _BEDTools: http://code.google.com/p/bedtools/
