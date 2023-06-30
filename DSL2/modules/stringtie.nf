process stringtie {
    tag "$name"
    label 'mega_memory'
    publishDir "${params.outdir}/star_mapped", pattern: "*{.gtf}*", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    container params.stringtie_container

    input:
    val(name)
    path(bam)
    path(bai)
    path(bw)
    each file(gtf)

    output:
    val(name)
//    file "${name}.gtf"
    file "${name}_for_DGE.gtf"


    script: 
    rf = params.stranded ? params.stranded == 'first-strand' ? '--rf' : '--fr' : ''
    """
    stringtie $bam -G $gtf -o ${name}.gtf $rf -a 8 -p $task.cpus
    stringtie $bam -G $gtf -o ${name}_for_DGE.gtf $rf -a 8 -e -p $task.cpus
    """
  }
