from optparse import OptionParser
from collections import Counter
from sets import Set
import os

def outputTracking(tracklines, outfn):
    print "Writing tracking file."
    outfile = open(outfn, 'w')
    written = 0
    samples = [str(x) + "_samp" for x in range(0, tracklines[0]['nsamples'])]
    for line in tracklines:
        samplist = [line[x] for x in samples]
        outline = "\t".join([line['transid'], line['locusid'],
                             line['ref_geneid'], line['ref_transid']] +
                            samplist) + "\n"
        outfile.write(outline)
        written = written + 1
    outfile.close()
    print "Wrote %d lines to %s." %(written, outfn)

def parseTracklineToDict(line):
    line = line.strip().split("\t")
    # first 4 columns are mandatory, the rest contain the sample info
    nsamples = len(line) - 4
    line = [nsamples] + line
    samples = [str(x) + "_samp" for x in range(0, nsamples)]
    keys = ["nsamples", "transid", "locusid", "ref_geneid",
            "ref_transid"] + samples
    
    linedict = dict(zip(keys, line))
    return linedict

def outputGTF(gtflines, outfn):
    print "Writing GTF file."
    outfile = open(outfn, 'w')
    written = 0
    for line in gtflines:
        
        outline = "\t".join([line['seqname'], line['source'],
                             line['feature'], line['start'],
                             line['end'], line['score'],
                             line['strand'], line['unknown'],
                             line['attribute']])
        outfile.write(outline)
        written = written + 1

    outfile.close()
    print "Wrote %d lines to %s." %(written, outfn)
    
def parseGTFlineToDict(line):
    values = line.split("\t")
    keys = ["seqname", "source", "feature", "start", "end", "score",
            "strand", "unknown", "attribute"]
    linedict = dict(zip(keys, values))
    return linedict

def parseTophatGTFlineToDict(line):
    linedict = parseGTFlineToDict(line)
    attributes = linedict["attribute"].split(";")
    del attributes[-1]
    attributes = [x.strip() for x in attributes]
    
    for attrib in attributes:
        attr = attrib.strip().split(" ")
        linedict[attr[0]] = attr[1].strip("\"")
    return linedict

def parseTophatGTF(gtf):
    print "Parsing GTF file."
    gtflines = []
    gfile = open(gtf, 'r')
    numlines = 0
    for line in gfile:
        numlines = numlines + 1
        gtflines.append(parseTophatGTFlineToDict(line))

    print "Processed %d lines in %s." %(numlines, gtf)
    gfile.close()

    return gtflines
    
def parseTracking(track):
    print "Parsing tracking file."
    tracklines = []
    tfile = open(track, 'r')
    numlines = 0
    for line in tfile:
        numlines = numlines + 1
        tracklines.append(parseTracklineToDict(line))

    print "Processed %d lines in %s." %(numlines, track)
    tfile.close()

    return tracklines
         
def filterGTFByTracking(gtflines, tracklines):
    idset = Set()
    kept = 0
    processed = 0
    for trackline in tracklines:
        idset.add(trackline['transid'])        

    newgtf = []
    for line in gtflines:
        processed = processed + 1
        if line['transcript_id'] in idset:
            newgtf.append(line)
            kept = kept + 1
        elif line['class_code'] == "=":
            newgtf.append(line)
            kept = kept + 1
            
    print "Kept %d out of %d lines in the GTF file." %(kept, processed)
    return(newgtf)

def filterGTFByExonCount(gtflines, threshold):
    print "Filtering for transcripts with at least %0.2f proportion of " \
          "exons in the longest transcript for a gene." %(threshold)
    transid_kept = Set([])
    transid_removed = Set([])
    genedict = {}
    total = 0
    kept = 0
    transid_removed = Set([])
    kept_lines = []
    
    # one trans_id per exon per line, so count up the number of
    # trans_ids to figure out the amount of exons attached to a gene_id
    for line in gtflines:
        total = total + 1
        genedict.setdefault(line['gene_id'], []).append(
            line['transcript_id'])

    print "Processed %d lines of the GTF file." %(total)
    
    # calculate the number of exons in each trans_id and save that id
    # if it passes the threshold
    for k,v in genedict.iteritems():
        genedict[k] = Counter(v)
        maxexons = float(genedict[k].most_common(1)[0][1])
        for trans_id, count in genedict[k].iteritems():
            if (count / maxexons) >= threshold:
                transid_kept.add(trans_id)
                
    for line in gtflines:
        if line['transcript_id'] in transid_kept:
            kept = kept + 1
            kept_lines.append(line)
        elif line['class_code'] == "=":
            kept = kept + 1
            transid_kept.add(line['transcript_id'])
            kept_lines.append(line)
        else:
            transid_removed.add(line['transcript_id'])

    print "Kept %d out of %d lines." %(kept, total)
    print "Kept %d transcripts." %(len(transid_kept))
    print "Removed %d transcripts." %(len(transid_removed))

    return kept_lines

def filterByNumsamples(tracklines, threshold):
    
    print "Filtering for lines with at least %d samples of " \
          "support." %(threshold)
    total = 0
    keptlines = []
    kept = 0
    for linedict in tracklines:
        total = total + 1
        if enoughSamples(linedict, threshold):
            kept = kept + 1
            keptlines.append(linedict)

    print "Kept %d lines out of %d lines had at least %d " \
          " samples of support " %(total, kept, threshold)
    return keptlines

def enoughSamples(linedict, threshold):
    """
    returns true if the number of samples supporting the transcript is
    greater than or equal to the threshold and false otherwise
    """

    support = linedict['nsamples'] - linedict.values().count("-")
    return(support >= threshold)
        
               
def main():
    parser = OptionParser()
    parser.add_option("-f", dest="gtf", default=None,
                      help="combined.gtf file from Cuffcompare")
    parser.add_option("-t", dest="tracking", default=None,
                      help="tracking file from Cuffcompare")
    parser.add_option("-s", dest="samples", default=0,
                      type="int",
                      help="number of samples a ID must appear in")
    parser.add_option("-e", dest="exon_thresh", default=0,
                      type="float",
                      help="minimum proportion of exons of longest isoform " \
                      "that must be in a transcript")

    (options, args) = parser.parse_args()

    if options.samples > 0:
        if options.tracking is None:
            print parser.print_help()
            exit(-1)
        trackout = options.tracking + ".filtered"
        if os.path.isfile(trackout):
            print "%s already exists, aborting." %(trackout)
            exit(-1)

        if options.gtf is None:
            print parser.print_help()
            exit(-1)
        gtfout = options.gtf + ".filtered"
        if os.path.isfile(gtfout):
            print "%s already exists, aborting." %(trackout)
            exit(-1)
            
    if options.exon_thresh > 0:
        if options.gtf is None:
            print parser.print_help()
            exit(-1)

        gtfout = options.gtf + ".filtered"
        if os.path.isfile(gtfout):
            print "%s already exists, aborting." %(gtfout)
            exit(-1)
            
    if options.samples + options.exon_thresh == 0.0:
        print parser.print_help()
        exit(-1)
        
    gtflines = parseTophatGTF(options.gtf)
    if options.samples > 0:
        tracklines = parseTracking(options.tracking)
        tracklines = filterByNumsamples(tracklines, int(options.samples))
        gtflines = filterGTFByTracking(gtflines, tracklines)
        outputTracking(tracklines, trackout)

    if options.exon_thresh > 0:
        gtflines = filterGTFByExonCount(gtflines, options.exon_thresh)


    outputGTF(gtflines, gtfout)
    
if __name__ == "__main__":
    main()
