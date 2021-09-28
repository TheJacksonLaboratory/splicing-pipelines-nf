# splicing-pipelines-nf
Anczukow-lab repository for the Nextflow splicing pipeline with `rMATS 4.1.0` (tag = v1.0) and `rMATS 4.1.1` (master).

## Introduction

The workflow processes raw data from FastQ inputs 
([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
    [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), aligns the reads
        ([STAR](https://github.com/alexdobin/STAR)), performs transcript assembly and quantification
            ([StringTie](https://ccb.jhu.edu/software/stringtie/)), detects alternative splicing events
                ([rMATs](http://rnaseq-mats.sourceforge.net/)), and generates an interactive QC report
                    ([MultiQC](http://multiqc.info/)).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. For example, it can be run on HPC (using [JAX](https://www.jax.org/)'s Sumner HPC) or on the cloud over ([Lifebit's CloudOS](https://lifebit.ai/cloudos) platform with AWS & GCloud). It comes with docker containers making installation trivial and results highly reproducible.

The input reads to this pipeline can come from 3 input sources:
![input_reads_graphic](https://raw.githubusercontent.com/lifebit-ai/images/master/jax_splicing/input_reads_graphic.png)

## Quick start on Sumner (JAX's HPC)

### Download latest changes to pipeline 

1) Go into splicing pipeline folder 
```bash
cd /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf
```

2) `status` to check branch and if there are unmerged changes. *Note*: `master` is the main branch. 

```bash
git status
```

3) `pull` any new changes. *Note*: you will need your github username and password

```bash
git pull
```

### Run test example to ensure the pipeline is working properly. 

#### Human test

1) Create new directory to run test and `cd` into that folder. Make this folder in individual's folder within `anczukow-lab`

2) Run test: 

```bash
sbatch /projects/anczukow-lab/splicing_pipeline/splicing-pipelines-nf/examples/human_test/human_test_main.pbs
```
3) To check progress: 

```bash
squeue -u [username]
```

or by looking at tail of output file: 

```bash
tail splicing.[jobnumber].out
```

### To see the help message
```
nextflow run main.nf --help
```

## Documentation

Documentation about the pipeline, found in the [`docs/`](docs) directory:

1. Intros to
    * [Containers](docs/containers.md)
    * [GitHub](docs/github.md)
2. [Pipeline overview](docs/pipeline_overview.md)
3. [Running the pipeline](docs/usage.md)
    * [Running on Sumner](docs/run_on_sumner.md)
    * [Running on CloudOS](docs/run_on_cloudos.md)
    * [Running locally](docs/run_locally.md)

## Pipeline DAG
<img src="https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/c23ee9552eb033dec087fb3b6fb01fe26716ce29/DAG.png" alt="splicing_pip_dag" align = "center" width="600"/>


## Changelog

### v 1.1 - Pipeline improvements
#### Fixes:
 - Adds missing trimmomatic logs to the multiqc report
 - Implemented correct support for input strandness in star process when `--stranded` is `second-strand` (was hardcoded to `strType=2` and only supported `first-strand` or `false` before)
#### Updates:
 - Updates the following tools:
   - **STAR** `2.7.3` -> `2.7.9a` NOTE: Requires a new index! (updated in test profile)
   - **Samtools** `1.10` -> `1.13`
   - **StringTie** `2.1.3b` -> `2.1.7`
   - **Gffread** `0.11.7` -> `0.12.7`
   - multiqc `1.8` -> `1.11`
   - deeptools `3.4.0` -> `3.5.1`
   - bioconductor-rtracklayer `1.46.0` -> `1.52.0`
   - gffcompare `0.11.2` -> `0.12.6`
   - bedtools `2.29.2` -> `2.30.0`
   - sra-tools `2.10.8` -> `2.11.0`
   - pigz `2.3.4` -> `2.6.0`
   - gdc-client `1.5.0` -> `1.6.1`
 - Moves all containers to https://hub.docker.com/u/anczukowlab

#### Maintenance:
 - Consideably reduces number of basic redundant CI tests by removing completely the `max_retries` matrix and `push` from `on: [push, pull_request]`
 - Adds CI test for sra-downloading pipeline pathway (only supported with docker profile for now)
#### QOL:
 - Adds saving of all the process .command* log files to results/process-logs folder
 - Adds pipeline workdir `--cleanup` option to clear all intermediate files on pipeline successful completion
 - Adds pipeline `--error_strategy` parameter to be able to specify pipeline error strategy directly from command line (doesn't work if specified in config linked by `-c` or `-config` nextflow params)

### v 1.0 - Initial pipeline release
