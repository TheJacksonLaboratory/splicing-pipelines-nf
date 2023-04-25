#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process star {
    tag "$name"
    label 'mega_memory'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'
    publishDir "${params.outdir}/star_mapped/${name}", pattern: "*{out.bam,out.bam.bai,out,ReadsPerGene.out.tab,SJ.out.tab}*" , mode: 'copy'
    publishDir "${params.outdir}/star_mapped/${name}", pattern: "*Unmapped*", mode: 'copy'
    publishDir "${params.outdir}/star_mapped/", mode: 'copy',
      saveAs: {filename -> 
          if (filename.indexOf(".bw") > 0) "all_bigwig/${name}.bw"
      }
    
    container params.star_container

    input:
    tuple val(name), path(reads), val(singleEnd)
    each path(index)
    each path(gtf)

    output:
    tuple val(name), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai") 
    path "*.out" 
    path "*ReadsPerGene.out.tab"
    path "*SJ.out.tab"
    path "*Log.out"
    path "*Unmapped*" optional true
    path "${name}.bw"

    script:
    out_filter_intron_motifs = params.stranded ? '' : '--outFilterIntronMotifs RemoveNoncanonicalUnannotated'
    out_sam_strand_field = params.stranded ? '' : '--outSAMstrandField intronMotif'
    xs_tag_cmd = params.stranded ? "samtools view -h ${name}.Aligned.sortedByCoord.out.bam | gawk -v q=${params.strType[params.stranded].strType} -f /usr/local/bin/tagXSstrandedData.awk | samtools view -bS - > Aligned.XS.bam && mv Aligned.XS.bam ${name}.Aligned.sortedByCoord.out.bam" : ''
    endsType = params.soft_clipping ? 'Local' : 'EndToEnd'
    star_mem = params.star_memory ? params.star_memory : task.memory
    avail_mem_bam_sort = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 2000000000}" : ''
    save_unmapped_reads = params.save_unmapped ? '--outReadsUnmapped Fastx' : ''

    """
    echo STAR \\
      --genomeDir ${index.toString().minus('.tar.gz')} \\
      --readFilesIn $reads \\
      --readMatesLengthsIn NotEqual \\
      --outFileNamePrefix ${name}. \\
      --runThreadN $task.cpus \\
      --readFilesCommand zcat \\
      --sjdbGTFfile $params.gtf \\
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
      --quantMode GeneCounts \\
      --outWigType None $out_filter_intron_motifs $out_sam_strand_field
      > ${name}.Aligned.sortedByCoord.out.bam
     """

workflow {

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
read_pairs_ch_singleEnd = Channel.from(params.singleEnd)
read_pairs_ch_index = Channel.fromPath(params.index)
read_pairs_ch_gtf = Channel.fromPath(params.gtf)

read_pairs_ch_star = read_pairs_ch.combine(read_pairs_ch_singleEnd)
read_pairs_ch_star.view()

star(read_pairs_ch_star, read_pairs_ch_index, read_pairs_ch_gtf)
}

}
