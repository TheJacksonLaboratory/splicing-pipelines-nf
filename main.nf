#!/usr/bin/env nextflow
/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'splicing-pipelines-nf' a pipeline repository to run Olga Anczukow's splicing pipeline.
 *
 * @authors
 * Marina Yurieva <marina.yurieva@jax.org>
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 * Carolyn Paisie
 * Phil Palmer <phil@lifebit.ai>
 * Olga Anczukow
 */

log.info "Splicing-pipelines - N F  ~  version 0.1"
log.info "====================================="
log.info "Reads                 : ${params.reads}"
log.info "Single-end            : ${params.singleEnd}"
log.info "Outdir                : ${params.outdir}"
log.info "\n"

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run jacksonlabs/splicing-pipelines-nf --reads '*_R{1,2}.fastq.gz' -profile docker
    Mandatory arguments:
       --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: docker, test and more.

    Generic:
      --singleEnd                   Specifies that the input is single-end reads
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

if (params.readPaths) {
  Channel
    .from(params.readPaths)
    .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
    .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
}
Channel
  .fromPath(params.adapter)
  .ifEmpty { exit 1, "Cannot find adapter sequence for trimming : ${params.adapter}" }
  .set { adapter }

/*--------------------------------------------------
  Run FastQC for quality control of input reads
---------------------------------------------------*/

process fastqc {
  tag "$name"
  publishDir "${params.outdir}/QC_raw", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_fastqc

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc --casava --threads $task.cpus $reads
  """
}

/*--------------------------------------------------
  Run Trimmomatic to trim input reads
---------------------------------------------------*/

process trimmomatic {
  tag "$name"
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_trimmomatic
  each file(adapter) from adapter

  output:
  file "*.fastq" into trimmed_reads

  script:
  mode = params.singleEnd ? 'SE' : 'PE'
  adapter_flag = params.adapter.endsWith("no_adapter.txt") ? '' : "ILLUMINACLIP:${adapter}:2:30:10:8:true"
  out = params.singleEnd ? "${name}_trimmed.fastq" : "${name}_paired_R1.fastq ${name}_unpaired_R1.fastq ${name}_paired_R2.fastq ${name}_unpaired_R2.fastq"
  """
  trimmomatic \
    $mode \
    -threads $task.cpus \
    -phred33 \
    $reads \
    $out \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:${params.readlength} \
    CROP:${params.readlength} $adapter_flag
  """
}