process prep_de {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/count_matrix", pattern: "{sample_lst.txt,*gene_count_matrix.csv,*transcript_count_matrix.csv}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    container params.stringtie_container  

    input:
    val(name)
    file(prepde_gtf)

    output:
    val(name)
    file "sample_lst.txt"
    file "*gene_count_matrix.csv"
    file "*transcript_count_matrix.csv"
    file("command-logs-*") optional true

    script:
    run_name = params.run_name ? params.run_name + "_" : ""
    date = new Date().format("MM-dd-yy")
    run_prefix = run_name + date

    """
    echo "${prepde_gtf.join("\n").toString().replace("_for_DGE.gtf", "")}" > samples.txt
    echo "${prepde_gtf.join("\n")}" > gtfs.txt
    paste -d ' ' samples.txt gtfs.txt > sample_lst.txt
    prepDE.py -i sample_lst.txt  -l $params.readlength \
              -g ${run_prefix}_gene_count_matrix.csv -t ${run_prefix}_transcript_count_matrix.csv
  
    ${params.savescript}
    """
   }
