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
