## Test the pipeline

See [here](../README.md#quick-start-on-sumner-jax-hpc-execution)

## Running your own analysis

The `main.pbs` executes an analysis comparing `MYC_MCF10A_0h` with 3 replicates and `MYC_MCF10_8h`.
The details of what needs to be configured to do this comparison analysis are found in three files:

All the analyses are kept in `analyses` subdirectory (by convention).   Encoded in the file name is the metadata details that outline the comparison that is being completed.  In this case capturing the statement above (`MYC_MCF10A_0h` vs `MYC_MCF10A_8h`).

These files can be specified via command line or via a config file.

To specify via command line:

* `--reads examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/reads.csv`
    This file contains the `sample_id` a short name uniquely defines the sample within this comparsion
    comma seperated with the complete path for the `left` and `right` `fastqs`.   
    
* `--b1 examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/b1.txt`
    This is a comma separated file containing 1 to many replicates for the `case` in the example.
    
* `--b2 examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/b2.txt`
    This is a comma separated file containing 1 to many replicates for the `control` in the example.