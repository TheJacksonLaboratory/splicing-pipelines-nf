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
