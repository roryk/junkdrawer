from sets import Set
import logging

def filterByMinLength(gtflines, size):
    newlines = []
    lengths = calculateLengths(gtflines)
    total = 0
    kept = 0
    for line in gtflines:
        total = total + 1
        if 'transcript_id' not in line:
            kept = kept + 1
            newlines.append(line)
            continue
        if line['transcript_id'] not in lengths:
            kept = kept + 1
            newlines.append(line)
            continue
        if lengths[line['transcript_id']] < size:
            continue
        kept = kept + 1
        newlines.append(line)

    logging.info("%d out of %d lines had a size greater than " \
                 "%d." %(kept, total, size))
    return newlines

def filterByMaxLength(gtflines, size):
    newlines = []
    lengths = calculateLengths(gtflines)
    total = 0
    kept = 0
    for line in gtflines:
        total = total + 1
        if 'transcript_id' not in line:
            kept = kept + 1
            newlines.append(line)
            continue
        if line['transcript_id'] not in lengths:
            kept = kept + 1
            newlines.append(line)
            continue
        if lengths[line['transcript_id']] > size:
            continue
        kept = kept + 1
        newlines.append(line)

    logging.info("%d out of %d lines had a size less than " \
                 "%d." %(kept, total, size))
    return newlines

def calculateLengths(gtflines):
    """
    calculate the lengths of each transcript in the gtf file by
    summing up the length of all the exons in the transcript
    """
    lengths = {}
    total_transcripts = 0
    for line in gtflines:
        if line['feature'] != "exon":
            continue
        size = abs(int(line['end']) - int(line['start'])) + 1
        if 'transcript_id' not in line:
            continue
        if line['transcript_id'] not in lengths:
            total_transcripts = total_transcripts + 1
            lengths[line['transcript_id']] = size
        else:
            lengths[line['transcript_id']] = lengths[line['transcript_id']] + \
                                             size

    logging.info("Processed %d transcripts." %(total_transcripts))
    return lengths

def outputLengths(lengths, outfn):
    outfile = open(outfn, 'w')
    written = 0
    for l in lengths.iteritems():
        written = written + 1
        l = [str(x) for x in l]
        outfile.write("\t".join(l) + "\n")
    logging.info("Wrote %d lines to %s." %(written, outfn))
    outfile.close()

def filterAttributes(gtflines, ffn):
    ffile = open(ffn, 'r')
    header = ffile.readline().strip()
    filter_set = Set()
    for line in ffile:
        filter_set.add(line.strip())

    newlines = []
    for line in gtflines:
        attrdict = attributeToDict(line)
        if not (header in attrdict):
            newlines.append(line)
            continue
        id = attrdict[header].replace("\"", "")
        if id in filter_set:
            newlines.append(line)
    
    return newlines

def reorderAttributes(gtflines, order):
    newlines = []
    order = order.split(" ")
    for line in gtflines:
        attrdict = attributeToDict(line)
        line['attribute'] = buildAttributeFieldFromDictWithOrder(attrdict,
                                                                 order)
        newlines.append(line)
    return(newlines)

def addAttribute(gtflines, afn):
    afile = open(afn, 'r')
    header = afile.readline()
    header = header.split("\t")
    logging.info("Attaching %s to %s." %(header[1], header[0]))
    linedict = {}
    for line in afile:
        line = line.split("\t")
        linedict[line[0]] = line[1]

    newlines = []
    for line in gtflines:
        attrdict = attributeToDict(line)
        if not (header[0] in attrdict):
            newlines.append(line)
            continue
        
        id = attrdict[header[0]].replace("\"", "")
        attrdict[header[1].strip()] = "\"" + linedict[id].strip() + "\""
        line['attribute'] = buildAttributeFieldFromDict(attrdict)
        newlines.append(line)

    return newlines

def swapAttributes(gtflines, source, replace):
    newlines = []
    repdict = dict(zip(replace, source))
    for line in gtflines:
        attrdict = attributeToDict(line)
        newdict = attrdict.copy()
        for r, s in repdict.items():
            if s in attrdict:
                newdict[r] = attrdict[s]
                line[r] = attrdict[s]
        line['attribute'] = buildAttributeFieldFromDict(newdict)
        newlines.append(line)
            
    return newlines

def delAttributes(gtflines, delete):
    newlines = []
    for line in gtflines:
        attrdict = attributeToDict(line)
        for x in delete:
            if x in attrdict:
                del attrdict[x]
            if x in line:
                del line[x]
            line['attribute'] = buildAttributeFieldFromDict(attrdict)
        newlines.append(line)
    return newlines

def GTFtoDict(infn):
    logging.info("Parsing the GTF file %s." %(infn))
    gtflines = []
    infile = open(infn, 'r')
    numlines = 0
    for line in infile:
        numlines = numlines + 1
        gtflines.append(parseGTFlineToDict(line))

    logging.info("Processed %d lines in %s." %(numlines, infn))
    infile.close()

    return gtflines

def buildAttributeFieldFromDict(attrdict):
    z = [x[0] + " " + x[1] for x in attrdict.items()]
    return("; ".join(z) + "\n")

def buildAttributeFieldFromDictWithOrder(attrdict, order):
    f = [0 for i in order]
    e = []
    for item in attrdict.items():
        if item[0] in order:
            f[order.index(item[0])] = item
        else:
            e.append(item)
    f = [x for x in f if x != 0]
    g = f + e
    z = [x[0] + " " + x[1] for x in g]
    return("; ".join(z) + "\n")
    
        
def attributeToDict(linedict):
    attributes = linedict["attribute"].split(";")
    attributes = filter(lambda x: x != "\n", attributes)

    return dict([x.strip().split(" ") for x in attributes])

def addAttributesToGTFline(linedict):
    attributes = linedict["attribute"].split(";")
    del attributes[-1]
    attributes = [x.strip() for x in attributes]
    
    for attrib in attributes:
        attr = attrib.strip().split(" ")
        linedict[attr[0]] = attr[1].strip("\"")
    return linedict

def parseGTFlineToDict(line):
    values = line.split("\t")
    keys = ["seqname", "source", "feature", "start", "end", "score",
            "strand", "unknown", "attribute"]
    linedict = dict(zip(keys, values))
    linedict = addAttributesToGTFline(linedict)
    
    return linedict

def outputGTF(gtflines, outfn):
    logging.info("Writing GTF file.")
    written = 0
    outfile = open(outfn, 'w')
    for line in gtflines:
        outline = formatGTFLine(line)
        outfile.write(outline)
        written = written + 1
        
    outfile.close()
    logging.info("Wrote %d lines to %s." %(written, outfn))

def formatGTFLine(line):
    outline = "\t".join([line['seqname'], line['source'],
                         line['feature'], line['start'],
                         line['end'], line['score'],
                         line['strand'], line['unknown'],
                         line['attribute']])
    return outline


def outputGTFout(gtflines):
    import sys
    logging.info("Writing GTF file.")
    written = 0
    for line in gtflines:
        outline = formatGTFLine(line)
        sys.stdout.write(outline)
        written = written + 1
    logging.info("Wrote %d lines." %(written))

def aggregateFeaturesByTranscript(gtflines):
    """
    aggregates a set of features by transcript_id
    if transcript_id does not exist it deletes the line
    sorts the features under a transcript by start site
    adds a start-end tuple to each transcript indicating where
    in the transcript each feature is (referenced to the mRNA, not
    the genome)
    """
    transcripts = {}
    for line in gtflines:
        if 'transcript_id' not in line:
            continue
        transcript_id = line['transcript_id']
        if transcript_id not in transcripts:
            transcripts[transcript_id] = [line]
        else:
            index = len(transcripts[transcript_id]) - 1

            while (transcripts[transcript_id][index]['start'] >
                   line['start']):
                index = index - 1
                if index == -1:
                    break
                
            transcripts[transcript_id].insert(index + 1, line)
    
    return transcripts

def aggregateFeaturesByGene(gtflines):
    genes = {}
    for line in gtflines:
        if 'gene_id' not in line:
            continue
        gene_id = line['gene_id']
        if gene_id not in genes:
            genes[gene_id] = [line]
        else:
            index = len(genes[gene_id]) - 1

            while (genes[gene_id][index]['start'] >
                   line['start']):
                index = index - 1
                if index == -1:
                    break

            genes[gene_id].insert(index + 1, line)

    return genes

def mergeOverlappedExons(genes):
    for gene in genes:
        exons = []
        for feature in genes[gene]:
            if feature['feature'] != "exon":
                continue
            if exons == []:
                exons.append(feature)
                continue
            
            for exon in reversed(exons):
                # the new feature is the furthest along, so add it and go to
                # the next feature
                if (exon['end'] < feature['start']):
                    exons.append(feature)
                    break
                # the old feature is contained in the new feature
                if ((exon['start'] >= feature['start']) and
                    exon['end'] <= feature['end']):
                    exon['end'] = feature['end']
                    exon['start'] = feature['start']
                    
                # the new feature is contained in the old feature
                # skip it
                if ((exon['start'] <= feature['start']) and
                    (exon['end'] >= feature['end'])):
                    break
                # the new feature overlaps with the old feature
                # merge them together
                if ((exon['start'] <= feature['start']) and
                    (exon['end'] <= feature['end'])):
                    exon['end'] = feature['end']
                    break
                if ((exon['end'] >= feature['end']) and
                    (exon['start'] >= feature['start'])):
                    exon['start'] = feature['start']
                    break
        genes[gene] = exons

    return genes

def addFeatureCoordinatesToTranscripts(transcripts):

    for trans_id in transcripts:
        start = 0
        for feature in transcripts[trans_id]:
            length = int(feature['end']) - int(feature['start']) + 1
            feature['fstart'] = start + 1
            feature['fend'] = start + length
            start = feature['fend']
            
    return transcripts

def orderFeaturesByTranscript(gtflines):
    """
    orders a set of features in each transcript by start site
    transcripts are not ordered
    """
    transcripts = aggregateFeaturesByTranscript(gtflines)
    new_gtflines = []
    for transcript in transcripts:
        for feature in transcript:
            new_gtflines.append(feature)
    return new_gtflines

def orderTranscriptsByChromosome(transcripts):
    """
    orders a set of transcripts by chromosome
    """
    chromosomes = {}
    for trans_id in transcripts.keys():
        seqname = transcripts[trans_id][0]['seqname']
        if seqname not in chromosomes:
            chromosomes[seqname] = [trans_id]
        else:
            index = len(chromosomes[seqname]) - 1
            ss = transcripts[chromosomes[seqname][index]][0]['start']

            while(int(ss) > int(transcripts[trans_id][0]['start'])):
                index = index - 1
                if index == -1:
                    break
                ss = transcripts[chromosomes[seqname][index]][0]['start']
            chromosomes[seqname].insert(index + 1, trans_id)

    return chromosomes

def unwindChromosomes(chromosomes, transcripts):
    gtflines = []

    chromosomelist = chromosomes.keys()
    chromosomelist.sort()

    for chromosome in chromosomelist:
        transcript_list = chromosomes[chromosome]
        for transcript_id in transcript_list:
            for gtfline in transcripts[transcript_id]:
                gtflines.append(gtfline)

    return gtflines
            


