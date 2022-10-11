#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import csv

def __main__():
    
    manifest = sys.argv[1]
    print("Input strandedness file:", manifest)

    count_2 = 0
    count_3 = 0
    count_4 = 0

    with open(manifest) as manifest_fh:
        data = csv.reader(manifest_fh, delimiter='\t')
        for line in data:
            count_2 += int(line[1])
            count_3 += int(line[2])
            count_4 += int(line[3])

    if count_4 > 200:
        print("Unstranded Data type")
        with open("infer_strandedness.txt", 'w') as fh:
            fh.write("false")
    else:
        if count_3 > count_4:
            print("first-strand")
            with open("infer_strandedness.txt", 'w') as fh:
                fh.write("first-strand")
        else:
            print("second-strand")
            with open("infer_strandedness.txt", 'w') as fh:
                fh.write("second-strand")


if __name__=="__main__": __main__()
