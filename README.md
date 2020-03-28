# splicing-pipelines-nf
Repository for the Nextflow splicing pipeline

## NextFlow Config Order

For the splicing pipeline we have:

### nextflow.config

* `nextflow.config` is the main config file found at the `root directory`. It includes defaults for all parameters, containers & links to other configs/profiles
* `base.config` to specify resources
* `igenomes.config` specifying the path to the igenomes resources available
    To mirror the setup on the CloudOS Universal Research Platform, adopting this standard setup.
    * Currently in `/projects/adeslat/igenomes/` GRCh38 has been downloaded
    * `STARIndex` downloaded from the AWS S3 Budget Amazon has donated to host these standards
    
* `MYC_MCF10A_0h_vs_MYC_MCF10A_8h.config`
    This file contains the configuration for the `Anczukow-lab` samples comparing MYC_MCF10A_0h triplicate samples with the MYC_MCF10A_8h triplicate samples.  As a convention, all the analyses are kept in `analyses` subdirectory.  Each config for comparison should contain 3 parameters
        * `reads` a `csv` file containing for paired read analysis
           * `sample_id`
           * `fastq1` left read
           * `fastq2` right read
        * `b1` the `txt` file a comma separated file containing 1 to many replicates used by the `rmats` program
        * `b2` the `txt` file a comma separated file containing 1 to many replicates used by the `rmats` program
    * examples directory
    * executors directory
        * sumner.conf - specifies the setup for HPC sumner
        *
igenomes.config to specify path to files for a given reference genome
test.config for using the test profile
user-defined configs eg the one you made for a specific analysis. This can be used to store a set of parameters and will be applied over the ones above if you use the -c option
finally there's command line arguments

## 
Download the pipeline and test it on a minimal dataset with a single command:
```bash
nextflow run jacksonlabs/splicing-pipelines-nf -profile test,<docker/sumner>
```

## Sumner Execution


```bash
nextflow run main.nf --reads examples/splicing_pipeline_testdata/single_end_reads.csv --genome GRCh38 -profile <docker/sumner>
```

For more info see how to [run pipelines on the Sumner HPC](docs/run_on_sumner.md)

## All parameters
```
Main arguments:
  --reads                       Path to input data CSV file specifying the reads sample_id and path to FASTQ files
  --singleEnd                   Specifies that the input is single-end reads
  --genome                      Name of iGenomes reference
  -profile                      Configuration profile to use. Can use multiple (comma separated)
                                Available: docker, test and more.

References:
  --gtf                         Path to GTF file
  --star_index                  Path to STAR index

Reads:
  --stranded                    Specifies that the input is stranded
  --b1                          Path to rMATs b1 file containing sample names
  --b2                          Path to rMATs b2 file containing sample names
  --adapter                     Path to adapter file
  --readlength                  Read length (default = 48)
  --overhang                    Overhang (default = readlength - 1)
  --mismatch                    Mismatch (default = 2)

Other:
  --max_cpus                    Maximum number of CPUs
  --max_memory                  Maximum memory
  --max_time                    Maximum time
  --skiprMATS                   Skip rMATS
  --skipMultiQC                 Skip MultiQC
  --outdir                      The output directory where the results will be saved
  --igenomes_base               Path to iGenomes base directory (default = 's3://ngi-igenomes/igenomes/')
```
