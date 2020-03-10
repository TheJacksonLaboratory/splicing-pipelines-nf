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
    .into { raw_reads_fastqc }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads_fastqc }
}

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