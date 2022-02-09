#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import pandas as pd

def __main__():
    
    manifest = sys.argv[1]
    reads = sys.argv[2]
    print("Input manifest file:", manifest)
    print("Input read file: ", reads)

    manifest_df = pd.read_csv(manifest, index_col=None, header=0, delimiter=",")

    if reads != "PASS":
        # process metadata
        reads_df = pd.read_csv(reads, index_col=None, header=0, delimiter=",")
        manifest_df = manifest_df[manifest_df['file_name'].isin(reads_df['file_name'].tolist())]

        if manifest_df.empty:
            print("Manifest file is empty after filtering.")
            sys.exit(404, "Manifest file is empty after filtering.")
        else:
            print("Number of samples in filtered manifest:")
            print(len(manifest_df))

    # save final manifest file
    manifest_df.to_csv("filtered_manifest.csv", sep=",", index=False) 

if __name__=="__main__": __main__()