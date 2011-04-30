# extract sequences for splice junctions mapped with tophat from the
# .juncs file
from optparse import OptionParser
from subprocess import call
import time
import os

def parseExtendedJuncLineToDict(line):
    linedict = {}
    bl = line.split('\t')
    idline = bl[4].split(':')
    linedict['chrom'] = idline[0]
    linedict['left'] = idline[1]
    linedict['right'] = idline[2]
    linedict['strand'] = idline[3]
    linedict['score'] = bl[5]
    return linedict

def parseTophatBedLineToDict(line):
    linedict = {}
    bl = line.split('\t')
    blocksizes = bl[10].split(',')
    blockstarts = bl[11].split(',')
    linedict['blockstarts'] = [int(x) for x in blockstarts]
    linedict['blocksizes'] = [int(x) for x in blocksizes]
    linedict['chromStart'] = int(bl[1])
    linedict['chromEnd'] = int(bl[2])
    linedict['name'] = bl[3]
    linedict['score'] = bl[4]
    linedict['strand'] = bl[5]
    linedict['chrom'] = bl[0]
    return(linedict)

def calculateJunctionEnds(linedict):
    """
    calculates the last base of the left end of the splice junction and the
    first base of the right end of the splice junction
    arguments: linedict a line dictionary
    returns: linedict augmented with 'left' and 'right' entries
    """
    left = linedict['chromStart'] + linedict['blocksizes'][0]
    right = linedict['chromStart'] + linedict['blockstarts'][1]

    linedict['left'] = left
    linedict['right'] = right
    return(linedict)

def formatJuncs(linedict, offset=0, ejuncs=False):
    """
    outputs a line in the juncs format
    """
    if not ejuncs:
        linedict = calculateJunctionEnds(linedict)
    chrom = linedict['chrom']
    left = str(int(linedict['left']) - offset)
    right = str(int(linedict['right']) + offset)
    strand = linedict['strand']
    score = str(linedict['score'])

    id = ":".join([chrom, left, right, strand])
    outline = "\t".join([chrom, left, right, strand, id, score])

    return outline

def formatExtendedJunc(linedict, offset):
    """
    returns two lines, one for each side of the junction
    """
    linedict = calculateJunctionEnds(linedict)
    chrom = linedict['chrom']
    left = str(linedict['left'])
    right = str(linedict['right'])
    left_off = str(linedict['left'] - offset)
    right_off = str(linedict['right'] + offset)
    strand = linedict['strand']
    score = str(linedict['score'])

    id = ":".join([chrom, left, right, strand])
    outleft = "\t".join([chrom, left_off, left, strand, id, score])
    outright = "\t".join([chrom, right, right_off, strand, id, score])
    return (outleft, outright)

def convertTophatBedToExtendedJuncs(infn, outfn, offset):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0

    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            linedict = parseTophatBedLineToDict(line)
            (outleft, outright) = formatExtendedJunc(linedict, offset)
            outfile.write(outleft + '\n')
            outfile.write(outright + '\n')
            
    infile.close()
    outfile.close()
    print "Processed %d junctions." %(junctions)

def convertTophatBedToJuncs(infn, outfn):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            linedict = parseTophatBedLineToDict(line)
            outline = formatJuncs(linedict)
            outfile.write(outline)

    infile.close()
    outfile.close()
    print "Processed %d junctions." %(junctions)    
            
    
def parseTophatJunctions(infn, outfn, keeplen):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    kept = 0
    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            bl = line.split('\t')
            blocksizes = bl[10].split(',')
            blockstarts = bl[11].split(',')
            blockstarts = [int(x) for x in blockstarts]
            blocksizes = [int(x) for x in blocksizes]
            bl[1] = int(bl[1])
            block1start = bl[1] + blockstarts[0]
            block1end = bl[1] + blockstarts[0] + blocksizes[0]
            block2start = bl[1] + blockstarts[1]
            block2end = bl[1] + blockstarts[1] + blocksizes[1]
            # trim the reads if the option to do that has been set
            if keeplen:
                block1start = block1end - keeplen
                block2end = block2start + keeplen
            if (block1start > block1end or block2start > block2end):
                continue
            kept = kept + 1
            outline = (bl[0] + "\t" + str(block1start) + "\t" + str(block1end) +
                       "\t" + bl[3] + "\t" + bl[4] + "\t" + bl[5])
            outfile.write(outline + '\n')
            outline = (bl[0] + "\t" + str(block2start) + "\t" + str(block2end) +
                       "\t" + bl[3] + "\t" + bl[4] + "\t" + bl[5])
            outfile.write(outline + '\n')
    infile.close()
    outfile.close()
    print "Read in %d junctions." %(junctions)
    print "Kept %d junctions. (%f)" %(kept, float(kept)/float(junctions))

def combineExtendedJuncs(infn, outfn, offset):
    """
    combines reads with the same junction name to reconstruct the
    original .junctions format of the file
    """
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    iddict = {}
    for line in infile:
        x = line.split("\t")
        iddict[x[4]] = parseExtendedJuncLineToDict(line)

    for value in iddict.values():
        outline = formatJuncs(value, offset, ejuncs=True)
        outfile.write(outline)

    infile.close()
    outfile.close()
                
        

def combineJunctionEnds(tmpfilename, outprefix):
    """
    combines reads with the same junction name to put together the
    whole splice junction sequence
    """
    tmpfilein = open(tmpfilename, 'r')
    juncdict = {}
    seqlen = 0
    junctions = 0

    for line in tmpfilein:
        if ">" in line:
            jname = line.split(">", 2)[1]
        else:
            seqlen = seqlen + len(line.rstrip())
            if jname not in juncdict:
                junctions = junctions + 1
                juncdict[jname] = line.rstrip()
            else:
                juncdict[jname] = juncdict[jname] + line.rstrip()

    tmpfilein.close()
    outfile = open(outprefix + ".fa", 'w')
    for key in juncdict:
        header = ">" + key 
        seq = juncdict[key] 
        outfile.write(header)
        outfile.write(seq)
        outfile.write("\n")
    outfile.close()

    print "Wrote %d junctions of %d total bases." %(junctions, seqlen)

def main():
    #usage = "usage: python %prog [options]"
    parser = OptionParser()
    parser.add_option("-f", dest="filename", default=None,
                      help="tophat junctions bed file")
    parser.add_option("-o", dest="outfile", 
                       help="output prefix")
    parser.add_option("-g", dest="genomefile", 
                       help="genome filename (FASTA format)")
    parser.add_option("-t", dest="keeplen", default=0,
                       help="keep only this many bases per side of junction")
    parser.add_option("-j", dest="juncs", default=False, action="store_true",
                      help="convert Tophat BED file to juncs file")
    parser.add_option("-c", dest="combine", default=False, action="store_true",
                      help="combine BED6 to BED 12 file by ID")

    (options, args) = parser.parse_args()
    
    if options.filename is None:
        print parser.print_help()
        exit(-1)
    if not os.path.isfile(options.filename):
        print "%s cannot be found." %(options.filename)
        exit(-1)
    if options.outfile is None:
        print parser.print_help()
        exit(-1)
    if options.genomefile is None:
        if not (options.juncs or options.combine):
            print parser.print_help()
            exit(-1)
    if not (options.juncs or options.combine):
        if not os.path.isfile(options.genomefile):
            print "%s cannot be found." %(options.genomefile)
            exit(-1)

    if options.juncs:
        if options.keeplen is 0:
            convertTophatBedToJuncs(options.filename, options.outfile)

        else:
            convertTophatBedToExtendedJuncs(options.filename, options.outfile,
                                           int(options.keeplen))
        exit(1)

    if options.combine:
        combineExtendedJuncs(options.filename, options.outfile,
                             int(options.keeplen))
        exit(1)
        
    tmpfilejuncs = "tmpjunctionfile" + str(time.time())
    parseTophatJunctions(options.filename, tmpfilejuncs, int(options.keeplen))
    tmpfileseq = "tmpseqfile" + str(time.time())
                                           
    print "Looking up sequences in %s with fastaFromBed." %(options.genomefile)
    call("fastaFromBed -fi %s -bed %s -name -fo %s" %(options.genomefile, tmpfilejuncs, tmpfileseq), shell=True)

    print "Combining junction ends."
    combineJunctionEnds(tmpfileseq, options.outfile)

    print "Cleaning up temporary files."
    os.remove(tmpfilejuncs)
    os.remove(tmpfileseq)

if __name__ == "__main__":
    main()
