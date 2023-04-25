#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//FastQC for quality control of input reads

process fastqc {

    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/QC/raw", pattern: "*_fastqc.{zip,html}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    container params.fastqc_container

    input:
    tuple val(name), file(reads)

    output:
    path "*_fastqc.{zip,html}"
    path ("command-logs-*") optional true

    script:
    """
    fastqc --casava --threads $task.cpus $reads

    ${params.savescript}
    """
  }

workflow {

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
read_pairs_ch.view()

Channel.fromFilePairs( params.reads, flat: true )
       .set{read_pairs_flat_ch}
read_pairs_flat_ch.view()

    fastqc( read_pairs_ch )
}
