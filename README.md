# splicing-pipelines-nf
Anczukow-lab repository for the Nextflow splicing pipeline with `rMATS 4.1.0` (tag = v1.0) and `rMATS 4.1.2` (master).

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

View changelog at [changelog.md](https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/changelog.md)
