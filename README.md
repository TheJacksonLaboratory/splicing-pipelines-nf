# splicing-pipelines-nf
Anczukow-lab repository for the Nextflow splicing pipeline with `rMATS 4.0.2`.

## Introduction

The workflow processes raw data from FastQ inputs 
([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
    [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), aligns the reads
        ([STAR](https://github.com/alexdobin/STAR)), perfroms transcript assembly and quantification
            ([StringTie](https://ccb.jhu.edu/software/stringtie/)), detects alterntive splicing events
                ([rMATs](http://rnaseq-mats.sourceforge.net/)), and generates an interactive QC report
                    ([MultiQC](http://multiqc.info/)).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. For example, it can be run on HPC (using [JAX](https://www.jax.org/)'s Sumner HPC) or on the cloud over ([Lifebit's CloudOS](https://lifebit.ai/cloudos) platfrom with AWS & soon GCloud). It comes with docker containers making installation trivial and results highly reproducible.

## Quick start on Sumner (JAX HPC) execution

1) Check out the source code from `github`.

```bash
git clone https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git
cd splicing-pipelines-nf
```

2) Then to execute on sumner, a user does the following bash script with the `slurm` `sbatch` command:

```bash
sbatch main.pbs
```

3) Progress can be checked using:

```bash
squeue | grep [username]
```

or by tailing the output file:

```bash
tail -f splicing.[jobnumber].out
```

## Documentation

Documentation about the pipeline, found in the docs/ directory:

1. [Pipeline overview](docs/pipeline_overview.md)
2. Pipeline configuration
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. Running the pipeline
    * [Running on Sumner](docs/run_on_sumner.md)
    * [Running on CloudOS](docs/run_on_cloudos.md)

## Diving into the details of the Example execution

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
  --gencode_gtf                 Path to gencode file (annotation gtf)

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
