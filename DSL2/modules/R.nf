process R {
   label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/", pattern: "{gffcmp.annotated.corrected.gff}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode:'copy'

    container params.R_container

    input:
    path(annotated_gtf)

    output:
    path "gffcmp.annotated.corrected.gff"

    script:
    """
    Rscript '/projects/anczukow-lab/yuriem/splicing_pipeline/splicing-pipelines-nf/bin/correct_gene_names.R'
    """
  }
