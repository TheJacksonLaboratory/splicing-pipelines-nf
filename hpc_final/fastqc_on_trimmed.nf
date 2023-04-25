#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fastqc_trimmed {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/QC/trimmed", pattern: "*_fastqc.{zip,html}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    container params.fastqc_container

    input:
    tuple val(name), path(reads), val(singleEnd)

    output:
    path "*_fastqc.{zip,html}"
    path("command-logs-*") optional true

    script:
    """
    fastqc --casava --threads $task.cpus $reads

    ${params.savescript}
    """
  }

workflow {

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
read_pairs_ch_singleEnd = Channel.from(params.singleEnd)

read_pairs_ch_fastqc = read_pairs_ch.combine(read_pairs_ch_singleEnd)

fastqc_trimmed( read_pairs_ch_fastqc )
}
