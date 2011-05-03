from sets import Set
import logging

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
        
        outline = "\t".join([line['seqname'], line['source'],
                             line['feature'], line['start'],
                             line['end'], line['score'],
                             line['strand'], line['unknown'],
                             line['attribute']])
        outfile.write(outline)
        written = written + 1
        
    outfile.close()
    logging.info("Wrote %d lines to %s." %(written, outfn))
