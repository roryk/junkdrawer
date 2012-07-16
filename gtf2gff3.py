#!/usr/bin/env python
import sys

def main(args):
    with open(args[0]) as in_handle:
        for line in in_handle.readlines():
            sline = line.split("\t")
            group = sline[8].split(";")
            out = sline[0:7] + [group[0].split(" ")[2].replace('\"', "")]
            print "\t".join(out)

if __name__ == "__main__":
    main(sys.argv[1:])
