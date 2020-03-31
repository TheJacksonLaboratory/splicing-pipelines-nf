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

## Quick start on Sumner (JAX's HPC)

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
3. [Running the pipeline](docs/usage.md)
    * [Running on Sumner](docs/run_on_sumner.md)
    * [Running on CloudOS](docs/run_on_cloudos.md)