#!/usr/bin/env nextflow
/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'splicing-pipelines-nf' a pipeline repository to run Olga Anczukow's splicing pipeline.
 *
 * @authors
 * Marina Yurieva <marina.yurieva@jax.org>
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 * Carolyn Paisie
 * Phil Palmer <phil@lifebit.ai>
 * Olga Anczukow
 */

// TODO: update logging & help message when pipeline is finished
log.info "Splicing-pipelines - N F  ~  version 0.1"
log.info "====================================="
log.info "Reads                 : ${params.reads}"
log.info "Single-end            : ${params.singleEnd}"
log.info "Outdir                : ${params.outdir}"
log.info "\n"

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run jacksonlabs/splicing-pipelines-nf --reads '*_R{1,2}.fastq.gz' -profile docker
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: docker, test and more.

    Generic:
      --singleEnd                   Specifies that the input is single-end reads
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false

if (params.singleEnd) {
  Channel
    .fromPath(params.reads)
    .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
    .splitCsv(skip:1)
    .map { sample_id, fastq -> [sample_id, file(fastq)] }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
} 
if (!params.singleEnd) {
  Channel
    .fromPath(params.reads)
    .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
    .splitCsv(skip:1)
    .map { sample_id, fastq1, fastq2 -> [ sample_id, [file(fastq1),file(fastq2)] ] }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
}
Channel
  .fromPath(params.adapter)
  .ifEmpty { exit 1, "Cannot find adapter sequence for trimming: ${params.adapter}" }
  .set { adapter }
Channel
  .fromPath(params.gtf)
  .ifEmpty { exit 1, "Cannot find GTF file: ${params.gtf}" }
  .into { gtf_star ; gtf_stringtie; gtf_stringtie_merge }
Channel
  .fromPath(params.star_index)
  .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
  .set { star_index }

/*--------------------------------------------------
  FastQC for quality control of input reads
---------------------------------------------------*/

process fastqc {
  tag "$name"
  publishDir "${params.outdir}/QC_raw", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_fastqc

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc --casava --threads $task.cpus $reads
  """
}

/*--------------------------------------------------
  Trimmomatic to trim input reads
---------------------------------------------------*/

process trimmomatic {
  tag "$name"
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_trimmomatic
  each file(adapter) from adapter

  output:
  set val(name), file(output_filename) into trimmed_reads

  script:
  mode = params.singleEnd ? 'SE' : 'PE'
  adapter_flag = params.adapter.endsWith("no_adapter.txt") ? '' : "ILLUMINACLIP:${adapter}:2:30:10:8:true"
  out = params.singleEnd ? "${name}_trimmed.fastq.gz" : "${name}_trimmed_R1.fastq.gz ${name}_unpaired_R1.fastq.gz ${name}_trimmed_R2.fastq.gz ${name}_unpaired_R2.fastq.gz"
  output_filename = params.singleEnd ? "${name}_trimmed.fastq.gz" : "${name}_trimmed_R{1,2}.fastq.gz"
  """
  trimmomatic \
    $mode \
    -threads $task.cpus \
    -phred33 \
    $reads \
    $out \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:${params.readlength} \
    CROP:${params.readlength} $adapter_flag
  """
}

/*--------------------------------------------------
  STAR to align trimmed reads
---------------------------------------------------*/

process star {
  tag "$name"
  publishDir "${params.outdir}/star_mapped", mode: 'copy'

  input:
  set val(name), file(reads) from trimmed_reads
  each file(index) from star_index
  each file(gtf) from gtf_star

  output:
  set val(name), file("${name}Aligned.sortedByCoord.out.bam"), file("${name}Aligned.sortedByCoord.out.bam.bai") into indexed_bam
  file "*.out" into alignment_logs
  file "*SJ.out.tab"
  file "*Log.out" into star_log
  file "*Unmapped*" optional true

  script:
  // TODO: check when to use `--outWigType wiggle` - for paired-end stranded stranded only?
  // TODO: check if `bw` file is in publishDir
  // TODO: test pipeline with paired-end data to ensure STAR has only two FASTQs to map
  // TODO: find a better solution to needing to use `chmod`
  out_filter_intron_motifs = params.stranded ? '' : '--outFilterIntronMotifs RemoveNoncanonicalUnannotated'
  overhang = params.overhang ? params.overhang : params.readlength - 1
  """
  # Decompress STAR index if compressed
  if [[ $index == *.tar.gz ]]; then
    tar -xvzf $index
  fi

  STAR \
    --genomeDir ${index.toString().minus('.tar.gz')} \
    --readFilesIn $reads \
    --readMatesLengthsIn NotEqual \
    --outFileNamePrefix $name \
    --runThreadN $task.cpus \
    --readFilesCommand zcat \
    --sjdbGTFfile $gtf \
    --sjdbOverhang $overhang \
     --alignSJoverhangMin 8 $out_filter_intron_motifs \
    --outFilterMismatchNmax $params.mismatch \
    --outFilterMultimapNmax 20 \
    --alignMatesGapMax 1000000 \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 100000000000 \
    --outBAMsortingThreadN $task.cpus \
    --outFilterType BySJout \
    --twopassMode Basic \
    --alignEndsType EndToEnd \
    --outWigType wiggle

  chmod a+rw $name*
  samtools index ${name}Aligned.sortedByCoord.out.bam
  bamCoverage -b ${name}Aligned.sortedByCoord.out.bam -o ${name}_old.bw 
  """
}

/*--------------------------------------------------
  Stringtie for transcript assembly and quantification 
---------------------------------------------------*/

process stringtie {
  tag "$name"
  publishDir "${params.outdir}/star_mapped", mode: 'copy'

  input:
  set val(name), file(bam), file(bam_index) from indexed_bam
  each file(gtf) from gtf_stringtie

  output:
  file "${name}.gtf" into stringtie_gtf
  file "${name}_for_DGE.gtf" into stringtie_deg_gtf

  script: 
  rf = params.stranded ? '--rf' : ''
  """
  stringtie $bam -G $gtf -o ${name}.gtf $rf -a 8 -p $task.cpus
  stringtie $bam -G $gtf -o ${name}_for_DGE.gtf $rf -a 8 -e -p $task.cpus
  """
}

/*--------------------------------------------------
  Stringtie merge GTF files
---------------------------------------------------*/

process stringtie_merge {
  publishDir "${params.outdir}/star_mapped", mode: 'copy'

  input:
  file('*.gtf') from stringtie_gtf.collect()
  file(gtf) from gtf_stringtie_merge

  output:
  file "stringtie_merged.gtf" into merged_gtf

  script:
  """
  ls -1 *.gtf > assembly_gtf_list.txt
  stringtie --merge -G $gtf -o stringtie_merged.gtf assembly_gtf_list.txt -p $task.cpus
  """
}