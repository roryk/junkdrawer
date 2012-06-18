#!/usr/local/bin/python
# shuffles the sequences in a fasta file
import random
from Bio import SeqIO
import fileinput

for record in SeqIO.parse(fileinput.input(), "fasta"):
    print ">" + record.id + "_shuffled"
    print ''.join(random.sample(record.seq, len(record.seq)))
