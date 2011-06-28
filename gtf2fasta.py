#!/usr/local/bin/python
from optparse import OptionParser
from BCBio import GFF
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import subprocess
import os
import sys

def lookupSequences(files):
    gtf_file = open(files['gtf_file'])
    records = []
    for rec in GFF.parse(gtf_file):
        chrom = rec.id
        for feature in rec.features:
            if feature.sub_features == []:
                seq = lookup_sequence(files, feature, chrom)
                id = feature.qualifiers['transcript_id'][0]

            else:
                seq = Seq("", generic_dna)
                id = feature.id
                for subf in feature.sub_features:
                    seq = seq + lookup_sequence(files, subf, chrom)

            records.append(SeqRecord(seq, id=id))
    SeqIO.write(records, sys.stdout, "fasta")
            
def lookup_sequence(files, feature, chrom):
    """
    use samtools to look up the sequence
    """
    args = [files['samtools'], "faidx", files['seq_file'], str(chrom) +
            ":" + str(feature.location.start) + "-" +
            str(feature.location.end)]

    child = subprocess.Popen(args, stdout=subprocess.PIPE)
    seq = ""
    for line in child.stdout:
        if line.strip()[0] == ">":
            continue
        seq = seq + line.strip()
    seq = Seq(seq, generic_dna)

    if feature.strand is 1:
        return seq
    elif feature.strand is -1:
        return seq.reverse_complement()
    else:
        print "strand is not 1 or -1, something is wrong"
        exit(-1)

def which(program):
    import os
    
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
                
    return None
                                                                    

def main():
    usage = "usage: gtf2fasta seq_file gtf_file"
    parser = OptionParser()
    (options, args) = parser.parse_args()
    samtools = which("samtools")
    
    if samtools is None:
        print "samtools must executable, add it to your path or " \
              "download it from http://samtools.sourceforge.net/"
        exit(-1)

    files = {}
    files['samtools'] = samtools

    if len(args) != 2:
        print usage
        exit(-1)

    files['seq_file'] = args[0]
    files['gtf_file'] = args[1]

    if not os.path.exists(files['seq_file']):
        print "seq_file does not exist"
        print usage
        exit(-1)

    if not os.path.exists(files['gtf_file']):
        print "gtf_file does not exist"
        print usage
        exit(-1)

    lookupSequences(files)
    
if __name__ == "__main__":
    main()
