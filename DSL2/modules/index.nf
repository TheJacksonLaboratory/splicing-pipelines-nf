process index {

    tag "$name"
    label 'mega_memory'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'
    publishDir "${params.outdir}/star_mapped", pattern: "{*.bai}" , mode: 'copy'

    container params.samtools_container

    input:
    val(name)
    path(bam)

    output:
    val(name)
    path(bam)
    path("*.bai")
  
    script:
    """
      samtools index ${bam}
    """
}
