#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import csv

def __main__():
    
    manifest = sys.argv[1]
    print("Input strandedness file:", manifest)

    if os.stat(manifest).st_size == 0:
        with open("infer_strandedness.txt", 'w') as fh:
            fh.write("false")
    else:
        with open(manifest) as manifest_fh:
            data = csv.reader(manifest_fh, delimiter='\t')
            for line in data:
                if line[0] == 'IU' or line[0] == 'U':
                    with open("infer_strandedness.txt", 'w') as fh:
                        fh.write("false")
                elif line[0] == 'ISR' or line[0] == 'SR':
                    with open("infer_strandedness.txt", 'w') as fh:
                        fh.write("first-strand")
                elif line[0] == 'ISF' or line[0] == 'SF':
                    with open("infer_strandedness.txt", 'w') as fh:
                        fh.write("second-strand")
                else:
                    with open("infer_strandedness.txt", 'w') as fh:
                        fh.write("false")

if __name__=="__main__": __main__()
