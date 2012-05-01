#!/usr/local/bin/python
import argparse
import sys
import logging
from juncs2seq import which
import subprocess

logging.basicConfig(level=logging.INFO)

IE_HEADER = ["exon_id", "inc", "exc"]
JUNC_HEADER = ["seqname", "left", "right", "id", "score", "strand"]
SEQNAMES = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
             'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr3', 'chr4', 'chr5',
             'chr6', 'chr7', 'chr8', 'chr9', 'chrUn', 'chrX']

def tsv2Dict(header, line):
    return dict(zip(header, line.split("\t")))

def countExonReads(exon, bamfile, samtools):
    args = [samtools, "view", bamfile.name, exon]
    child = subprocess.Popen(args, stdout=subprocess.PIPE)
    reads = 0
    for line in child.stdout:
        reads = reads + 1

    return reads

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ie_file', nargs='?', type=argparse.FileType('r'),
                        help="a inclusion/exclusion file.")
    parser.add_argument('j1_file', nargs='?', type=argparse.FileType('r'),
                        help="a file in juncs format")
    parser.add_argument('bamfile', nargs='?', type=argparse.FileType('r'),
                       help="a file in bam format")
    args = parser.parse_args()

    samtools = which("samtools")

    for seq in SEQNAMES:

        logging.info("Reading in %s from event file." %(seq))
        # ditch the header
        args.ie_file.readline()
        events = []
        for line in args.ie_file:
            if (not line.count("#")) and (line.split(":")[0] == seq):
                events.append(IEEvent(line))
        logging.info("Read in %d events." %(len(events)))
                
        logging.info("Reading %s from junction file." %(seq))
        juncs = {}
        for line in args.j1_file:
            seqid = line.split("\t")[0]
            if seqid == seq:
                junc = Junction(line)
                juncs[junc.get_id()] = junc
        logging.info("Read in %d junctions on %s." %(len(juncs), seq))

        out_count = 0
        for event in events:
            # just to make sure
            if event.seqname == seq:
                out_count = out_count + 1
                if(out_count % 1000 == 0):
                    logging.info("Processed %d events" %(out_count))
                event.count_inclusion(juncs)
                event.count_exclusion(juncs)
                # add all of the reads that map to the exon
                event.inc_counts = event.inc_counts + countExonReads(event.exon_id, args.bamfile, samtools)
                event.output_counts()
        logging.info("Output %d event counts." %(out_count))

        # reset the two files
        args.j1_file.seek(0)
        args.ie_file.seek(0)
        
    args.j1_file.close()
    args.ie_file.close()

class IEEvent:
    def __init__(self, line):
        self.event = dict(zip(IE_HEADER, line.split("\t")))
        self.exon_id = self.event["exon_id"]
        self.inc_id = self.event["inc"].split(",")
        self.exc_id = self.event["exc"].split(",")
        self.inc_counts = 0
        self.exc_counts = 0
        self.seqname = line.split(":")[0]

    def count_inclusion(self, juncs):
        for inc in self.inc_id:
            if inc in juncs:
                self.inc_counts = self.inc_counts + juncs[inc].get_score()
        return self.inc_counts

    def count_exclusion(self, juncs):
        for exc in self.exc_id:
            if exc in juncs:
                self.exc_counts = self.exc_counts + juncs[exc].get_score()
        return self.exc_counts

    def output_counts(self):
        sys.stdout.write("\t".join([self.exon_id, str(self.inc_counts),
                                    str(self.exc_counts)]) + "\n")

class Junction:
    def __init__(self, line):
        self.junction = dict(zip(JUNC_HEADER, line.split("\t")))

    def get_id(self):
        return self.junction["id"]

    def get_score(self):
        return int(self.junction["score"])

    def __eq__(self, other):
        if isinstance(other, Junction):
            return self.get_id() == other.get_id()
        elif isinstance(other, str):
            return self.get_id() == other
        else:
            return False

if __name__ == "__main__":
    main()
