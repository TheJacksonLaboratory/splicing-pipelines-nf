process gffread {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/", pattern: "{gffcmp.annotated.corrected.gtf}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode:'copy'

    container params.gffread_container

    input:
    path(gffread_gtf)

    output:
    path "gffcmp.annotated.corrected.gtf"
    path("command-logs-*") optional true

    script:
    """
    gffread -E gffcmp.annotated.corrected.gff -T -o gffcmp.annotated.corrected.gtf
    ${params.savescript}
    """
  }
