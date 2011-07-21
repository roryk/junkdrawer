import re
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from gtfUtils import *
import os
import subprocess
import logging


def calculate_transcript_length(transcript):
    """
    using a list of features in a transcript, calculate the length
    of the transcript
    """
    max_length = 0
    for x in transcript:
        if 'fend' not in x:
            continue
        if x['fend'] > max_length:
            max_length = x['fend']

    return max_length
            

def add_start_and_stop_codons(orf, transcript):
    """
    add start and stop codons to a transcript
    will only work for eurkaryotic transcripts
    """
    
    # grab a template for the start and stop codons, model it on an
    # exon, but removing the exon_specific fields
    tmp = ""
    for x in transcript:
        if x['feature'] == "exon":
            tmp = x.copy()
            break
    if tmp == "":
        return
    tmp = delAttributes([tmp], ["exon_number"])[0]

    length = calculate_transcript_length(transcript)

    # on the minus strand, the transcript is going the other way, but
    # we had to reverse complement it to determine where the ORF is,
    # this flips the coordinates back so it points at the right place
    # on the - strand
    if tmp['strand'] == "-":
        orf['start'] = length - orf['start'] + 1 - 3
        orf['end'] = length - orf['end'] + 1 + 3
    # don't bother to annotate a transcript without a sane strand
    elif tmp['strand'] != "+":
        return transcript


    def make_codon(orf, transcript, feature_type, new_template):

        template = new_template.copy()
        # if the strand information isnt set, skip this transcript
        if (template['strand'] != "+") and (template['strand'] != "-"):
            return transcript
        
        if (feature_type == "start_codon"):
            offset = 3
                
            template['fstart'] = orf['start']
            template['fend'] = orf['start'] + offset
            overlapped = find_overlapping_features([orf['start']], transcript)
            if overlapped == []:
                return []
            overlapped = overlapped[0]
            template['start'] = int(overlapped['start']) + \
                                   int(orf['start']) - \
                                   int(overlapped['fstart']) + 1
            template['end'] = template['start'] + offset
                                                             
        elif feature_type == "stop_codon":
            offset = -3
                
            template['fend'] = orf['end']
            template['fstart'] = orf['end'] + offset
            overlapped = find_overlapping_features([orf['end']], transcript)
            if overlapped == []:
                return []
            overlapped = overlapped[0]
            template['end'] = int(overlapped['start']) + \
                                int(orf['end']) - \
                                int(overlapped['fstart']) + 1
            template['start'] = template['end'] + offset

        template['feature'] = feature_type
        template['source'] = "gtfutils"

        return template

    start_codon = make_codon(orf, transcript, "start_codon", tmp)
    stop_codon = make_codon(orf, transcript, "stop_codon", tmp)

    # if there is no start and stop codon, then its not likely a real open
    # reading frame, so don't annotate it with one and not the other
    if (start_codon == []) or (stop_codon == []): 
        return transcript

    transcript.append(start_codon)
    transcript.append(stop_codon)

    return transcript
    
def find_longest_orf(seq):
    """
    finds the longest open reading frame in a cDNA sequence string
    returns a dictionary with the key-value pairs:
    longest: length of longest ORF
    seq: cDNA sequence of ORF
    start: start codon position in the cDNA sequence
    end: stop codon position in the cDNA sequence
    """
    r = re.compile("(ATG(\w\w\w)*?(TAG|TAA|TGA))", re.IGNORECASE)
    z = r.finditer(seq)
    
    d = {'longest': 0, 'start': 0, 'end': 0, 'seq': ""}
    for item in z:
        if len(item.group()) > d['longest']:
            d['longest'] = len(item.group())
            d['start'] = item.start()
            d['end'] = item.end()
            d['seq'] = item.group()

    # set the frame to be 1 based not 0 based
    d['start'] = d['start'] + 1
    d['end'] = d['end']
    
    return d

def get_sequence_of_transcript(samtools, fasta_file, transcript):
    seq = Seq("", generic_dna)
    for feature in transcript:
        if feature['feature'] != "exon":
            continue
        feature_seq = get_sequence_of_feature(samtools, fasta_file, feature)
        seq = seq + Seq(feature_seq, generic_dna)
    strand = feature['strand']
    if strand == "-":
        seq = seq.reverse_complement()
    return seq.tostring()


def get_sequence_of_feature(samtools, fasta_file, feature):
    args = [samtools, "faidx", fasta_file,
            feature['seqname'] + ":" + str(feature['start']) +
            "-" + str(feature['end'])]
    child = subprocess.Popen(args, stdout=subprocess.PIPE)
    seq = ""
    for line in child.stdout:
        if line.strip()[0] == ">":
            continue
        seq = seq + line.strip()
    return seq

def find_overlapping_features(positions, transcript):
    """
    find the features that overlap a dictionary of features in a
    transcript in a list of positions
    """
    overlapped = []
    for feature in transcript:
        for position in positions:
            if (position >= feature['fstart']) and (position <= feature['fend']):
                overlapped.append(feature)
                break

    return overlapped

def which(program):
    
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
    logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',
                        level=logging.INFO)
    description = "Add estimated start and stop codons to a GTF file."
    parser = ArgumentParser(description = description)
    parser.add_argument("fasta_file", metavar="seq_file", type=str, 
                        help="sequence of transcripts in fasta format")
    parser.add_argument("gtf_file", metavar="gtf_file", type=str, 
                        help="gtf file")

    samtools = which("samtools")

    if samtools is None:
        print "samtools must be executable, add it to your path or " \
              "download it from http://samtools.sourceforge.net/"
        exit(-1)
    
    args = parser.parse_args()

    if not os.path.exists(args.gtf_file):
        parser.help()
        exit(-1)
    if not os.path.exists(args.fasta_file):
        parser.help()
        exit(-1)

    gtflines = GTFtoDict(args.gtf_file)
    
    x = aggregateFeaturesByTranscript(gtflines)
    x = addFeatureCoordinatesToTranscripts(x)
    for transcript in x:
        seq = get_sequence_of_transcript(samtools, args.fasta_file, x[transcript])
        orf = find_longest_orf(seq)
        z = add_start_and_stop_codons(orf, x[transcript])
        outputGTFout(z)
        
    #z = orderTranscriptsByChromosome(x)
    #gtflines = unwindChromosomes(z, x)
    #outputGTFout(gtflines)

if __name__ == "__main__":
    main()
