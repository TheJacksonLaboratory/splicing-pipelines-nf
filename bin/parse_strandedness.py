#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

def __main__():
    
    manifest = sys.argv[1]
    threshold = float(sys.argv[2])
    print("Input strandedness file:", manifest)

    fraction1 = 0
    fraction2 = 0

    with open(manifest) as manifest_fh:
        for line in manifest_fh:
            line = line.split(':')
            if line[0]=='Fraction of reads explained by "1++,1--,2+-,2-+"':
                fraction1 = float(line[1].strip())
            elif line[0]=='Fraction of reads explained by "1+-,1-+,2++,2--"':
                fraction2 = float(line[1].strip())
            else:
                pass

    if 1 - threshold < (fraction1 / fraction2) < 1 + threshold:
        print("Unstranded Data type")
        with open("infer_strandedness.txt", w) as fh:
            fh.write("false")
    else:
        if fraction1 > fraction2:
            print("first-strand")
            with open("infer_strandedness.txt", w) as fh:
                fh.write("first-strand")
        else:
            print("second-strand")
            with open("infer_strandedness.txt", w) as fh:
                fh.write("second-strand")


if __name__=="__main__": __main__()
