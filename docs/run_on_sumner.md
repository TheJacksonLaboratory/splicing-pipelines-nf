## Test the pipeline

See [here](../README.md##quick-start-on-sumner-jaxs-hpc)

## Running your own analysis

### 1. Create input files

1) Create `reads` input CSV file
    - You will need to create a CSV file containing the path to your input `reads`. You can see examples for [single-end](../examples/testdata/single_end/reads.csv) and [paired-end](../../examples/testdata/paired_end/paired_end_reads_replicates.csv) data
2) Optional: create the `b1` and `b2` text files used by rMATS.
    - As the pipeline takes FASTQ (not BAM) input the values will need to be the `sample_id` as specified in the [`reads`](../examples/testdata/paired_end/paired_end_reads_replicates.csv) file. See example [`b1.txt`](../examples/testdata/single_end/b1.txt) and [`b2.txt`](../examples/testdata/single_end/b2.txt)

### 2. Setup your own configuration file

Here you will be [adding your own custom config](https://nf-co.re/usage/configuration#custom-configuration-files)

The config file will likely be specific to your user and analysis. **You do not need to edit the pipeline code to configure the pipeline**.

To create your own custom config (to specify your input parameters) you can copy and edit this [`example.config`](../conf/examples/MYC_MCF10A_0h_vs_MYC_MCF10A_8h_sumner.config) file.

The file contains all available parameters with sensible defaults which you can find more information on [here](usage.md#all-available-parameters). You will need to specify the path to your `reads` and `b1`/`b2` input files. This string can be a relative path from the directory which you run Nextflow in, an absolute path or even a link as shown by the [`example.config`](../conf/examples/MYC_MCF10A_0h_vs_MYC_MCF10A_8h_sumner.config).

### 3. Run the pipeline

Once you have created the input files and config then you can run the Nextflow pipeline using the [`main.pbs`](../main.pbs) script. You will need to modify the `main.pbs` script so that it runs the pipeline with your config profile and the `sumner` profile, for example:
```
./nextflow run main.nf -config /path/to/my_file.config -profile sumner -resume
```

Then [run the `main.pbs` script](../README.md#quick-start-on-sumner-jaxs-hpc) to submit jobs to the cluster.

### Examples

See [`MYC_MCF10A_0h_vs_MYC_MCF10A_8h.config`](../conf/examples/MYC_MCF10A_0h_vs_MYC_MCF10A_8h_sumner.config) for an example analysis comparing the 0h and 8h timepoints

The analysis compares `MYC_MCF10A_0h` with 3 replicates and `MYC_MCF10_8h`.
The details of what needs to be configured to do this comparison analysis are found in three files:

All the analyses can be kept in `examples/analyses` subdirectory (by convention). Encoded in the file name is the metadata details that outline the comparison that is being completed.  In this case capturing the statement above (`MYC_MCF10A_0h` vs `MYC_MCF10A_8h`).

These files can be specified via command line or via a config file.

To specify via command line:

* `--reads examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/reads.csv`
    This file contains the `sample_id` a short name uniquely defines the sample within this comparsion
    comma seperated with the complete path for the `left` and `right` `fastqs`.   
    
* `--b1 examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/b1.txt`
    This is a comma separated file containing 1 to many replicates for the `case` in the example.
    
* `--b2 examples/analyses/MYC_MCF10A_0h_vs_MYC_MCF10A_8h/b2.txt`
    This is a comma separated file containing 1 to many replicates for the `control` in the example.