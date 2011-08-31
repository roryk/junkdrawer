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
    this is 0 based.
    arguments: linedict a line dictionary
    returns: linedict augmented with 'left' and 'right' entries
    """
    left = linedict['chromStart'] + linedict['blocksizes'][0]
    right = linedict['chromStart'] + linedict['blockstarts'][1]

    # -1 because it is 0 based
    linedict['left'] = left
    linedict['right'] = right + 1
    return(linedict)
    
