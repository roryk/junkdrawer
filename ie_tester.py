#!/usr/local/bin/python

import numpy as np
import scipy.stats as stats
import sys

def parseIELine(line):
    if line[0] == "#":
        return []
    else:
        x = line.split("\t")
        z = [int(y) for y in x[1:]]
        return [x[0]] + z

def main():
    header = sys.stdin.readline().split("\t")
    for line in sys.stdin:
        x = parseIELine(line)
        z = [(y < 5) for y in x[1:]]
        if True in z:
            continue
        z = np.array([[x[1], x[2]], [x[3], x[4]]])
        c = stats.chi2_contingency(z, correction=False)
        if c[1] < 0.01:
            print line
            print c

if __name__ == "__main__":
    main()
