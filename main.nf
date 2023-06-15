#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Check if user has set adapter sequence. If not set is based on the value of the singleEnd parameter
adapter_file = params.adapter ? params.adapter : params.singleEnd ? "$baseDir/adapters/TruSeq3-SE.fa" : "$baseDir/adapters/TruSeq3-PE.fa"

// Set overhang to read length -1, unless user specified
overhang = params.overhang ? params.overhang : params.readlength - 1

// if download_from was specified
download_from = params.download_from ? params.download_from : ""
key_file = params.key_file ? params.key_file : "$baseDir/examples/assets/no_key_file.txt"

// Get minlength, if user does not specify, set to read length
minlen = params.minlen ? params.minlen : params.readlength

// if minlength != read length, turn on variable read length 
variable_read_length = minlen == params.readlength ? false : true

// get run name and date prefix for counts matrix and multiqc
run_name = params.run_name ? params.run_name + "_" : ""
date = new Date().format("MM-dd-yy")
run_prefix = run_name + date

//FastQC for quality control of input reads

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

process trimmomatic {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/trim", pattern: "*_R{1,2}.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    container params.trimmomatic_container

    input:
    tuple val(name), path(reads), val(singleEnd), path(adapter)

    output:
    tuple val(name), path(output_filename), val(singleEnd)
    path ("logs/${name}_trimmomatic.log")
    path ("command-logs-*") optional true

    script:
    mode = singleEnd ? 'SE': 'PE'
    out = singleEnd ? "${name}_trimmed.fastq.gz" : "${name}_trimmed_R1.fastq.gz ${name}_unpaired_R1.fastq.gz ${name}_trimmed_R2.fastq.gz ${name}_unpaired_R2.fastq.gz"
    output_filename = singleEnd ? "${name}_trimmed.fastq.gz" : "${name}_trimmed_R{1,2}.fastq.gz"
    slidingwindow = params.slidingwindow ? 'SLIDINGWINDOW:4:15' : ''
    keepbothreads = singleEnd == true ? '' : ':2:true'
    
    """
    trimmomatic \
      $mode \
      -threads $task.cpus \
      -phred33 \
      $reads \
      $out \
      ILLUMINACLIP:${adapter}:2:30:10${keepbothreads} \
      LEADING:3 \
      TRAILING:3 \
      $slidingwindow \
      MINLEN:${params.minlen} \
      CROP:${params.readlength}
    mkdir logs
    cp .command.log logs/${name}_trimmomatic.log

    ${params.savescript}
    """
  }

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

if(params.debug) {
    pre_script_run_resource_status = """
        echo ==========================
        echo Debug Summary
        echo ==========================
        echo === Machine type ===
        uname --all 
        echo === Machine memory ===
        free -g -t
        echo === Machine CPU ===
        nproc --all
        echo === Pre-script run disk-free ===
        df -h /
        echo ==========================
        """
    post_script_run_resource_status = """
        echo === Post-script run disk-free ===
        df -h /
        echo ==========================
        """
  }else{
    pre_script_run_resource_status = ""
    post_script_run_resource_status = ""
  }

process star {
    tag "$name"
    label 'mega_memory'

    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'
    publishDir "${params.outdir}/star_mapped", pattern: "*{out.bam,out,ReadsPerGene.out.tab,SJ.out.tab}*" , mode: 'copy'    
    
    container params.star_container

    input:
    tuple val(name), path(reads)
    val(singleEnd)
    each path(index)
    each path(gtf)

    output:
    val(name)
    path("${name}.Aligned.sortedByCoord.out.bam")

    script:
    out_filter_intron_motifs = params.stranded ? '' : '--outFilterIntronMotifs RemoveNoncanonicalUnannotated'
    out_sam_strand_field = params.stranded ? '' : '--outSAMstrandField intronMotif'
    endsType = params.soft_clipping ? 'Local' : 'EndToEnd'
    star_mem = params.star_memory ? params.star_memory : task.memory
    avail_mem_bam_sort = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 2000000000}" : ''
    save_unmapped_reads = params.save_unmapped ? '--outReadsUnmapped Fastx' : ''
    
    """
    ${pre_script_run_resource_status}
    
     if [[ $index == *.tar.gz ]]; then
      tar -xvzf $index
    fi
    
   STAR \\
      --genomeDir ${index.toString().minus('.tar.gz')} \\
      --readFilesIn $reads \\
      --readMatesLengthsIn NotEqual \\
      --outFileNamePrefix ${name}. \\
      --runThreadN $task.cpus \\
      --readFilesCommand zcat \\
      --sjdbGTFfile $gtf \\
      --sjdbOverhang $params.overhang \\
      --alignSJDBoverhangMin $params.sjdbOverhangMin \\
      --outFilterScoreMinOverLread $params.filterScore \\
      --outFilterMatchNminOverLread $params.filterScore \\
      --outFilterMismatchNmax $params.mismatch \\
      --outFilterMultimapNmax 20 \\
      --alignMatesGapMax 1000000 \\
      --outSAMattributes All \\
      --outSAMtype BAM SortedByCoordinate \\
      $avail_mem_bam_sort \\
      --outBAMsortingThreadN $task.cpus \\
      --outFilterType BySJout \\
      --twopassMode Basic \\
     --alignEndsType $endsType \\
      --alignIntronMax 1000000 \\
      $save_unmapped_reads \\
      --quantMode GeneCounts \\
      --outWigType None $out_filter_intron_motifs $out_sam_strand_field
      """
  }

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
    """
    echo "${prepde_gtf.join("\n").toString().replace("_for_DGE.gtf", "")}" > samples.txt
    echo "${prepde_gtf.join("\n")}" > gtfs.txt
    paste -d ' ' samples.txt gtfs.txt > sample_lst.txt
    prepDE.py -i sample_lst.txt  -l $params.readlength \
              -g ${run_prefix}_gene_count_matrix.csv -t ${run_prefix}_transcript_count_matrix.csv
  
    ${params.savescript}
    """
   }

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

 
workflow {

reads = Channel.fromPath( params.reads)
  .splitCsv(skip:1)
  .map { sample_id, fastq1, fastq2 -> [ sample_id, [file(fastq1),file(fastq2)] ] }


if (params.bams) {
  Channel
    .fromPath(params.bams)
    .ifEmpty { exit 1, "Cannot find BAMs csv file : ${params.bams}" }
    .splitCsv(skip:1)
    .map { name, bam, bai -> [ name, file(bam), file(bai) ] }
//    .into { indexed_bam; indexed_bam_rmats }
    .view()
} 

fastqc( reads )

reads_trimmomatic = Channel.fromFilePairs( params.trimmomatic_reads, checkIfExists: true )
reads_trimmomatic_singleEnd = Channel.from(params.singleEnd)
reads_trimmomatic_adapter = Channel.fromPath(params.adapter)

reads_trimmomatic_trim = reads_trimmomatic.combine(reads_trimmomatic_singleEnd)
reads_trimmomatic_final = reads_trimmomatic_trim.combine(reads_trimmomatic_adapter)

trimmomatic(reads_trimmomatic_final)

reads_fastqc_trimmed = Channel.fromFilePairs( params.fastqc_on_trimmed_reads, checkIfExists: true )
reads_fastqc_trimmed_singleEnd = Channel.from(params.singleEnd)
reads_fastqc_trimmed_final = reads_fastqc_trimmed.combine(reads_fastqc_trimmed_singleEnd)

fastqc_trimmed( reads_fastqc_trimmed_final)

reads_star = Channel.fromFilePairs( params.star_reads, checkIfExists: true )
reads_star_singleEnd = Channel.from(params.singleEnd)
reads_star_index = Channel.fromPath(params.index)
reads_star_gtf = Channel.fromPath(params.gtf)
reads_final_star = reads_star.combine(reads_star_singleEnd)

(name, bam, bai, bw) = star(reads, reads_final_star, reads_star_index, reads_star_gtf) | index 
bamcov_result = bamCoverage(index.out)
stringtie_result = stringtie(bamCoverage.out, reads_star_gtf)
prepde_result = prep_de(stringtie.out)
stringtie_merge_result = stringtie_merge(stringtie.out, reads_star_gtf)

gffcompare_gtf = Channel.fromPath(params.stringtie_merge_gtf)
gffcompare_result = gffcompare(gffcompare_gtf, reads_star_gtf)
R_gtf = Channel.fromPath(params.annotated_gtf)
R_result = R(R_gtf) | gffread
}
