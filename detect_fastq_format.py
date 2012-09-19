#!/usr/bin/env python
import sys

_FASTQ_RANGES = {"sanger": [33, 73],
                 "solexa": [59, 104],
                 "illumina_1.3+": [64, 104],
                 "illumina_1.5+": [66, 104],
                 "illumina_1.8+": [33, 74]}


def detect_fastq_format(in_file, MAX_RECORDS=1000000):
    """
    detects the format of a fastq file
    will return multiple formats if it could be more than one
    """
    kept = list(_FASTQ_RANGES.keys())
    with open(in_file) as in_handle:
        records_read = 0
        for i, line in enumerate(in_handle):
            # get the quality line
            if records_read >= MAX_RECORDS:
                break
            if i % 4 is 3:
                records_read += 1
                for c in line:
                    formats = kept
                    if len(formats) == 1:
                        return formats
                    for form in formats:
                        if (_FASTQ_RANGES[form][0] > ord(c) or
                            _FASTQ_RANGES[form][1] < ord(c)):
                            kept.remove(form)

    return formats


if __name__ == "__main__":
    if len(sys.argv) is not 1:
        print "usage: detect_fastq_format.py fastq_file"
    print "Possible formats %s" % (str(detect_fastq_format(sys.argv[1])))
