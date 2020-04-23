#!/usr/bin/env nextflow
/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'splicing-pipelines-nf' a pipeline repository to run Olga Anczukow's splicing pipeline.
 *
 * @authors
 * Laura Urbanski <laura.urbanski@jax.org> first author of the Post-processing portion of the workflow!
 * Marina Yurieva <marina.yurieva@jax.org>
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 * Carolyn Paisie
 * Phil Palmer <phil@lifebit.ai>
 * Anne Deslattes Mays <adeslatt@gmail.com>
 * Olga Anczukow
 */

log.info "Splicing-pipelines - N F  ~  version 0.1"
log.info "====================================="
log.info "Assembly name         : ${params.assembly_name}"
log.info "Reads                 : ${params.reads}"
log.info "Single-end            : ${params.singleEnd}"
log.info "GTF                   : ${params.gtf}"
log.info "STAR index            : ${params.star_index}"
log.info "Stranded              : ${params.stranded}"
log.info "rMATS pairs file      : ${params.rmats_pairs ? params.rmats_pairs : 'Not provided'}"
log.info "Adapter               : ${params.adapter.endsWith('no_adapter.txt') ? 'Not provided' : params.adapter}"
log.info "Read Length           : ${params.readlength}"
log.info "Overhang              : ${params.overhang}"
log.info "Mismatch              : ${params.mismatch}"
log.info "Outdir                : ${params.outdir}"
log.info "Max CPUs              : ${params.max_cpus}"
log.info "Max memory            : ${params.max_memory}"
log.info "Max time              : ${params.max_time}"
log.info ""
log.info "\n"

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run jacksonlabs/splicing-pipelines-nf --reads my_reads.csv --gtf genome.gtf --star_index star_dir -profile base,sumner
    
    Main arguments:
      --reads                       Path to input data CSV file specifying the reads sample_id and path to FASTQ files
      --gtf                         Path to GTF file
      --star_index                  Path to STAR index
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: base, docker, sumner, test and more.

    Reads:
      --rmats_pairs                 Path to file containing b1 & b2 samples names space seperated, one row for each rMATS comparison
      --singleEnd                   Specifies that the input is single-end reads
      --stranded                    Specifies that the input is stranded
      --adapter                     Path to adapter file
      --readlength                  Read length (default = 48)
      --overhang                    Overhang (default = readlength - 1)
      --mismatch                    Mismatch (default = 2)

    Other:
      --assembly_name               Genome assembly name
      --max_cpus                    Maximum number of CPUs
      --max_memory                  Maximum memory
      --max_time                    Maximum time
      --skiprMATS                   Skip rMATS
      --skipMultiQC                 Skip MultiQC
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

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
// Loaded as a file rather than a channel incase params.adapter is undefined
adapter = file(params.adapter)
Channel
  .from(params.assembly_name)
  .ifEmpty { exit 1, "Genome assembly name not set"}
  .set { assembly_name }

Channel
  .fromPath(params.gtf)
  .ifEmpty { exit 1, "Cannot find GTF file: ${params.gtf}" }
  .into { gtf_star ; gtf_stringtie; gtf_stringtie_merge; gtf_rmats }

Channel
  .fromPath(params.star_index)
  .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
  .set { star_index }
if (params.rmats_pairs) {
  Channel
    .fromPath(params.rmats_pairs)
    .ifEmpty { exit 1, "Cannot find rMATS pairs file : ${params.rmats_pairs}" }
    .splitCsv(sep:' ')
    .map { row -> 
      def b1 = row[0].toString().split(',')
      def b2 = row[1].toString().split(',')
      [ b1 + b2]
    }
    .set { samples}
}

/*--------------------------------------------------
  FastQC for quality control of input reads
---------------------------------------------------*/

process fastqc {
  tag "$name"
  label 'process_medium'
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
  label 'low_memory'
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
  label 'high_memory'
  publishDir "${params.outdir}/star_mapped", mode: 'copy'

  input:
  set val(name), file(reads) from trimmed_reads
  each file(index) from star_index
  each file(gtf) from gtf_star

  output:
  set val(name), file("${name}Aligned.sortedByCoord.out.bam"), file("${name}Aligned.sortedByCoord.out.bam.bai") into (indexed_bam, indexed_bam_rmats)
  file "*.out" into alignment_logs
  file "*SJ.out.tab"
  file "*Log.out" into star_log
  file "*Unmapped*" optional true
  file "${name}_old.bw"

  script:
  // TODO: check when to use `--outWigType wiggle` - for paired-end stranded stranded only?
  // TODO: find a better solution to needing to use `chmod`
  out_filter_intron_motifs = params.stranded ? '' : '--outFilterIntronMotifs RemoveNoncanonicalUnannotated'
  out_sam_strand_field = params.stranded ? '' : '--outSAMstrandField intronMotif'
  overhang = params.overhang ? params.overhang : params.readlength - 1
  xs_tag_cmd = params.stranded ? "samtools view -h ${name}Aligned.sortedByCoord.out.bam | awk -v strType=2 -f /usr/local/bin/tagXSstrandedData.awk | samtools view -bS - > Aligned.XS.bam && mv Aligned.XS.bam ${name}Aligned.sortedByCoord.out.bam" : ''
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
    --alignSJoverhangMin 8 \
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
    --alignIntronMax 1000000 \
    --outReadsUnmapped Fastx \
    --outWigType wiggle $out_filter_intron_motifs $out_sam_strand_field

  chmod a+rw $name*
  $xs_tag_cmd
  samtools index ${name}Aligned.sortedByCoord.out.bam
  bamCoverage -b ${name}Aligned.sortedByCoord.out.bam -o ${name}_old.bw 
  """
}

/*--------------------------------------------------
  Stringtie for transcript assembly and quantification 
---------------------------------------------------*/

process stringtie {
  tag "$name"
  label 'process_medium'
  publishDir "${params.outdir}/star_mapped", mode: 'copy'

  input:
  set val(name), file(bam), file(bam_index) from indexed_bam
  each file(gtf) from gtf_stringtie

  output:
  file "${name}.gtf" into stringtie_gtf
  file "${name}_for_DGE.gtf" into stringtie_dge_gtf

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
  label 'process_medium'
  publishDir "${params.outdir}/star_mapped", mode: 'copy'

  input:
  file('*.gtf') from stringtie_gtf.collect()
  file(gtf) from gtf_stringtie_merge

  output:
  file "stringtie_merged.gtf" into merged_gtf
  file "gffcmp.*" into gffcmp

  script:
  """
  ls -1 *.gtf > assembly_gtf_list.txt
  stringtie --merge -G $gtf -o stringtie_merged.gtf assembly_gtf_list.txt -p $task.cpus
  gffcompare -R -V -r $gtf stringtie_merged.gtf
  correct_gene_names.R
  gffread -E gffcmp.annotated.corrected.gff -T -o gffcmp.annotated.corrected.gtf
  """
}

/*--------------------------------------------------
  rMATS to detect alternative splicing events
---------------------------------------------------*/

if (params.rmats_pairs) {

  indexed_bam_rmats
    .map { name, bam, bai -> bam }
    .set { bam }
  
  // Group BAMs for each rMATS execution
  samples
    .map { row -> 
      // Create unique row id by joining all the sample names
      def row_id = (row[0]).join(",").toString().replace(",", "")
      def samples_id = []
      row[0].each { sample ->
        samples_id.add([sample, row_id])
      }
      samples_id
    }
    .flatMap()
    .combine(bam, by:0)
    .map { sample_id, row_id, bam -> [row_id, bam] }
    .groupTuple()
    .map { row_id, bams -> bams }
    .set { bams }

  process rmats {
    label 'high_memory'
    publishDir "${params.outdir}/rMATS_out", mode: 'copy'

    when:
    !params.skiprMATS

    input:
    file(bams) from bams
    file(gtf) from gtf_rmats
    file (gffcmp) from gffcmp
    val(fasta) from assembly_name


    output:
    file "*"

    script:
    mode = params.singleEnd ? 'single' : 'paired'
    n_samples_replicates = bams.size()
    n_replicates = n_samples_replicates.intdiv(2)
    bam_groups = bams.collate(n_replicates)
    b1_bams = bam_groups[0].join(",")
    b2_bams = bam_groups[1].join(",")
    """
    echo $b1_bams > b1.txt
    echo $b2_bams > b2.txt
    rmats.py --b1 b1.txt --b2 b2.txt --gtf $gtf --od ./ -t $mode --nthread $task.cpus --readLength ${params.readlength}
    rmats_config="config_for_rmats_and_postprocessing.txt"
    echo b1 b1.txt > \$rmats_config
    echo b2 b2.txt >> \$rmats_config
    echo rmats_gtf $gtf >> \$rmats_config
    echo rmats_gtf       ${gtf} >> \$rmats_config
    echo ref_gtf         ${gtf} >> \$rmats_config
    echo fasta           ${params.fasta} >> \$rmats_config
    echo reads           ${params.singleEnd ? 'single' : 'paired'} >> \$rmats_config
    echo readlen         ${params.readlength} >> \$rmats_config
    
    LU_postprocessing.R
    """
  }

} else {

  indexed_bam_rmats
    .map { name, bam, bai -> [name, bam] }
    .toSortedList { entry -> entry[0] }
    .flatten()
    .collate(4, false)
    .set { paired_samples }

  process paired_rmats {
    tag "$name1 $name2"
    label 'high_memory'
    publishDir "${params.outdir}/rMATS_out", mode: 'copy'

    when:
    !params.skiprMATS

    input:
    set val(name1), file(bam1), val(name2), file(bam2) from paired_samples
    file (gtf) from gtf_rmats
    file (gffcmp) from gffcmp

    output:
    file "*"

    script:
    mode = params.singleEnd ? 'single' : 'paired'
    """
    ls $bam1 > b1.txt
    ls $bam2 > b2.txt
    rmats.py --b1 b1.txt --b2 b2.txt --gtf $gtf --od ./ -t $mode --nthread $task.cpus --readLength ${params.readlength}

    rmats_config="config_for_rmats_and_postprocessing.txt"
    echo b1 b1.txt > \$rmats_config
    echo b2 b2.txt >> \$rmats_config
    echo rmats_gtf $gtf >> \$rmats_config
    echo rmats_gtf       ${gtf} >> \$rmats_config
    echo ref_gtf         ${gtf} >> \$rmats_config
    echo fasta           ${params.fasta} >> \$rmats_config
    echo reads           ${params.singleEnd ? 'single' : 'paired'} >> \$rmats_config
    echo readlen         ${params.readlength} >> \$rmats_config
    
    LU_postprocessing.R
    """
  }
}

/*--------------------------------------------------
  MultiQC to generate a QC HTML report
---------------------------------------------------*/

process multiqc {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
  file ('alignment/*') from alignment_logs.collect().ifEmpty([])

  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data"

  script:
  """
  multiqc . -m fastqc -m star
  """
}