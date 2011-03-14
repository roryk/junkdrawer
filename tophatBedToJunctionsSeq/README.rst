======================
tophatBedToJunctionSeq
======================
Takes a Tophat .juncs file and a genome in FASTA format and outputs
the sequences making up the junction.

Requires that fastaFromBed be installed from [[BEDTools][http://code.google.com/p/bedtools/]


Usage: tophatBedToJunctionsSeq.py [options]

Options:
  -h, --help     show this help message and exit
  -f FILENAME    tophat junctions bed file
  -o OUTFILE     output prefix
  -g GENOMEFILE  genome filename (FASTA format)
