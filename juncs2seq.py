#!/usr/local/bin/python

import sys
import logging
import os
import subprocess
import string

logging.basicConfig(level=logging.INFO)

JUNC_HEADER = ["seqname", "left", "right", "id", "score", "strand"]

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

def _parse_junc(line):
    return dict(zip(JUNC_HEADER, line.split("\t")))

def _calculate_ends(junction, keep, ignore):
    if keep >= 0:
        left_start = int(junction['left']) - keep + 1
        left_end = int(junction['left']) + 1
        right_start = int(junction['right']) + 1
        right_end = int(junction['right']) + keep + 1

        if ignore:
            left_end = left_end - 1
            right_start = right_start + 1
        
    if keep < 0:
        left_start = int(junction['left']) + 1
        left_end = int(junction['left']) - keep + 1
        right_start = int(junction['right']) + keep + 1
        right_end = int(junction['right']) + 1
        
        if ignore:
            left_start = left_start + 1
            right_end = right_end - 1

    left = str(left_start) + "-" + str(left_end)
    right = str(right_start) + "-" + str(right_end)

    return (left, right)

def _build_query(junction, location):
    query = junction['seqname'] + ":" + location
    return query

def _reverse_complement(seq):
    table = string.maketrans("ACGT", "TGCA")
    seq = seq.translate(table)[::-1]
    return seq

def _lookup_junctions(files, queries):
    """
    looks up both sides of nseqs number of junctions and concatenates
    them together
    """
    args = [files['samtools'], "faidx", files['seq_file']] + queries
    child = subprocess.check_output(args).split("\n")
    
    seqs = []
    seq = ""
    begin = True
    # combine every other sequence together since the queries are in
    # pairs for junctions
    for line in child:
        if line == "":
            continue
        if line.strip()[0] == ">":
            continue
        if begin:
            seq = line.strip()
            begin = False
        else:
            seq = seq + line.strip()
            seqs.append(seq)
            begin = True

    return seqs

def read_juncs():
    juncs = []
    njuncs = 0
    for line in sys.stdin:
        njuncs = njuncs + 1
        juncs.append(_parse_junc(line))

    logging.info("Read %d junctions." %(njuncs))
    
    return juncs

def read_juncs_nice(CHUNK_SIZE):
    juncs = []
    i = 0
    for line in sys.stdin:
        i = i + 1
        if i > CHUNK_SIZE:
            juncs.append(_parse_junc(line))
            return juncs
        juncs.append(_parse_junc(line))
        
    return juncs

def _build_output_lines(out1, out2, rt):
    name = ">" + out1 + "\n"
    if rt and (out1.strip()[-1] == "-"):
        seq = _reverse_complement(out2) + "\n"
    else:
        seq = out2 + "\n"
    return name + seq

def main():
    
    description = "Converts a set of junction to sequences."
    usage = "juncs2seq genome N\n\t" \
            "genome: genome file to use (fasta format and indexed)\n\t" \
            "N: keep N nucleotides on each side of the junction. " \
            "This can be a negative number which will get a piece of the intervening intron.\n\t" \
            "--i: ignore the actual junction nucleotide\n\t" \
            "--rt: reverse complement sequences on the minus strand."
    
    if len(sys.argv) < 3:
        print description
        print usage
        exit(1)

    CHUNK_SIZE = 1000
    files = {'seq_file': sys.argv[1], 'samtools': which('samtools')}
    keep = int(sys.argv[2])

    ignore = False
    if sys.argv.count("--i"):
        ignore = True

    rt = False
    if sys.argv.count("--rt"):
        rt = True

    queries = []
    names = []
    seqs = []
    written = 0
    logging.info("Processing junctions.")
    juncs = read_juncs_nice(CHUNK_SIZE)

    while juncs != []:
        for junc in juncs:
            names.append(junc['id'])
            loc = _calculate_ends(junc, keep, ignore)
            queries.append(_build_query(junc, loc[0]))
            queries.append(_build_query(junc, loc[1]))
        
        seqs = seqs + _lookup_junctions(files, queries)
        outputs = zip(names, seqs)
        # reset back to the initial state
        queries = []
        seqs = []
        names = []
        for output in outputs:
            sys.stdout.write(_build_output_lines(output[0], output[1], rt))
            written = written + 1
        # output a progress indicator since this can take a while
        sys.stderr.write(".")    
        juncs = read_juncs_nice(CHUNK_SIZE)
    sys.stderr.write("\n")
    
    logging.info("Wrote %d junctions total." %(written))

if __name__ == "__main__":
    main()
