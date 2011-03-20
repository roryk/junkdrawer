======================
tophatBedToJunctionSeq
======================

Takes a Tophat .juncs file and a genome in FASTA format and outputs
the sequences making up the junction.

Requires that fastaFromBed from BEDTools_ be installed and executable.

Usage: tophatBedToJunctionsSeq.py [options]

Options:
  -h, --help     show this help message and exit
  -f FILENAME    tophat junctions bed file
  -o OUTFILE     output prefix
  -g GENOMEFILE  genome filename (FASTA format)
  -t SIZE	 # of bases to keep on each side of splice junction
     		 this is optional
.. _BEDTools: http://code.google.com/p/bedtools/
