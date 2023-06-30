process bamCoverage {

    tag "$name"
    label 'mega_memory'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'
    publishDir "${params.outdir}/star_mapped", pattern: "{*.bw}", mode: 'copy'

    container params.deeptools_container

    input:
    val(name)
    path(bam)
    path(bai)


    output:
    val(name)
    path(bam)
    path(bai)
    path ("*.bw") , emit: bigwig, optional: true

    script:
  
    """
    bamCoverage -b ${bam} -o ${bai}.bw
    """
}
