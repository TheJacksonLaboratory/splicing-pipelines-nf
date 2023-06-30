process stringtie_merge {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped", pattern: "*{.gtf,.txt}*", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    container params.stringtie_container

    input:
    val(name)
    path(prepde_gtf)
    path(gtf)

    output:
    val(name)
    path "stringtie_merged.gtf"
    path "assembly_gtf_list.txt"
    path("command-logs-*") optional true

    script:
    """
    ls -1 *.gtf > assembly_gtf_list.txt
    stringtie --merge -G $gtf -o stringtie_merged.gtf assembly_gtf_list.txt -p $task.cpus
    ${params.savescript}
    """
  }
