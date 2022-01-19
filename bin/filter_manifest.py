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
            sys.exit(404, "Manifest file is empty after filtering.")
        
        if len(reads_df[~reads_df['file_name'].isin([manifest_df['file_name'].tolist()])])>0:
            print("The following file_name IDs where not found in manifest:")
            print(reads_df[~reads_df['file_name'].isin([manifest_df['file_name'].tolist()])])
            reads_df[~reads_df['file_name'].isin([manifest_df['file_name'].tolist()])].to_csv("not_found_GTEX_samples.txt", index=False)

    # save final manifest file
    manifest_df.to_csv("filtered_manifest.csv", sep=",", index=False) 

if __name__=="__main__": __main__()