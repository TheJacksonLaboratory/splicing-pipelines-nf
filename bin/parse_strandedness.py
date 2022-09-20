#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

def __main__():
    
    manifest = sys.argv[1]
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

    if 1 - 0.2 < (fraction1 / fraction2) < 1.2:
        print("Unstranded Data type")
    else:
        if fraction1 > fraction2:
            print("first-strand")
        else:
            print("second-strand")


if __name__=="__main__": __main__()
