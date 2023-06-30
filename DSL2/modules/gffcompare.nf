process gffcompare {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/", pattern: "{gffcmp.*}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode:'copy'

    container params.gffcompare_container

    input:
    path(gtf)
    path(stringtie_merge_gtf)

    output:
    path "gffcmp.*"
    path("command-logs-*") optional true

    script:
    """
    gffcompare -R -V -r $gtf stringtie_merged.gtf
    ${params.savescript}
    """
 }
