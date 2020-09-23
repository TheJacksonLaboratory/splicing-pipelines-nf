#!/usr/bin/env nextflow
/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'splicing-pipelines-nf' a pipeline repository to run Olga Anczukow's splicing pipeline.
 *
 * @authors
 * Laura Urbanski <laura.urbanski@jax.org> first author of the Post-processing portion of the workflow!
 * Marina Yurieva <marina.yurieva@jax.org>
 * Brittany Angarola <brittany.angarola@jax.org>
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 * Carolyn Paisie
 * Phil Palmer <phil@lifebit.ai>
 * Sangram Sahu <sangram@lifebit.ai>
 * Anne Deslattes Mays <adeslatt@gmail.com>
 * Olga Anczukow
 */

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --reads my_reads.csv --gtf genome.gtf --star_index star_dir -profile base,sumner
    
    Input files:
      --reads                       Path to reads.csv file, which specifies the sample_id and path to FASTQ files for each read or read pair (path).
                                    This file is used if starting at beginning of pipeline. 
                                    (default: no reads.csv)
      --bams                        Path to bams.csv file which specifies sample_id and path to BAM and BAM.bai files (path)
                                    This file is used if starting pipeline at Stringtie.
                                    (default: no bams.csv)
      --rmats_pairs                 Path to rmats_pairs.txt file containing b1 (and b2) samples names (path)
                                    (default: no rmats_pairs specified) 
      --run_name                    User specified name used as prefix for output files
                                    (defaut: no prefix, only date and time)
      --download_from               Database to download FASTQ/BAMs from (available = 'TCGA', 'GTEX' or 'SRA') (string)
                                    (default: false)
      --key_file                    For downloading reads, use TCGA authentication token (TCGA) or dbGAP repository key (GTEx, path)
                                    (default: false)  
 
    Main arguments:
      --gtf                         Path to reference GTF file (path)
                                    (default: no gtf specified) 
      --assembly_name               Genome assembly name (available = 'GRCh38' or 'GRCm38', string)
                                    (default: false)
      --star_index                  Path to STAR index (path)
                                    (default: read length)
      --singleEnd                   Specifies that the input is single-end reads (bool)
                                    (default: false)
      --stranded                    Specifies that the input is stranded ('first-strand', 'second-strand', false (aka unstranded))
                                    (default: 'first-strand')
      --readlength                  Read length - Note that all reads will be cropped to this length(int)
                                   (default: no read length specified)
      -profile                      Configuration profile to use. Can use multiple (comma separated, string)
                                    Available: base, docker, sumner, test and more.

    Trimmomatic: 
      --minlen                      Drop the read if it is below a specified length (int)
                                    Default parameters turn on --variable-readlength
                                    To crop all reads and turn off, set minlen = readlength (NOTE: this will turn off soft clipping)                                
                                    (default: 20)
      --slidingwindow               Perform a sliding window trimming approach (bool)
                                    (default: true)
      --adapter                     Path to adapter file (path)  
                                    (default: TruSeq3 for either PE or SE, see singleEnd parameter)
                                
    Star:                    
      --mismatch                    Number of allowed mismatches per read (SE) or combined read (PE) (int)
                                    SE ex. read length of 50, allow 2 mismatches per 50 bp
                                    PE ex. read length of 50, allow 2 mismatches per 100 bp 
                                    (default: 2)
      --overhang                    Overhang (int)
                                    (default: readlength - 1)
      --filterScore                 Controls --outFilterScoreMinOverLread and outFilterMatchNminOverLread
                                    (default: 0.66)
      --sjdbOverhangMin              Controls --alignSJDBoverhangMin (int)
                                    (default: 3)
      --star_memory                 Max memory to be used by STAR to sort BAM files.
                                    (default: Available task memory)

    rMATS:                              
      --statoff                     Skip the statistical analysis (bool)
                                    If using only b1 as input, this must be turned on.
                                    (default: false)
      --paired_stats                Use the paired stats model (bool)
                                    (default: false)
      --novelSS                     Enable detection of unnanotated splice sites (bool)
                                    (default: false)
      --mil                         Minimum Intron Length. Only impacts --novelSS behavior (int)
                                    (default: 50)
      --mel                         Maximum Exon Length. Only impacts --novelSS behavior (int)
                                    (default: 500)

    Other:
      --test                        For running trim test (bool)
                                    (default: false)
      --max_cpus                    Maximum number of CPUs (int)
                                    (default: ?)  
      --max_memory                  Maximum memory (memory unit)
                                    (default: 80)
      --max_time                    Maximum time (time unit)
                                    (default: ?)
      --skiprMATS                   Skip rMATS (bool)
                                    (default: false)
      --skipMultiQC                 Skip MultiQC (bool)
                                    (default: false)
      --outdir                      The output directory where the results will be saved (string)
                                    (default: directory where you submit the job)
      --mega_time                   Sets time limit for processes withLabel 'mega_memory' in the main.nf using the base.config (time unit)     
                                    (default: 20.h)
      --gc_disk_size                Only specific to google-cloud executor. Adds disk-space for few aggregative processes.
                                    (default: "200 GB" based on 100 samples. Simply add 2 x Number of Samples)

    See here for more info: https://github.com/TheJacksonLaboratory/splicing-pipelines-nf/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Check if read length is set
if (!params.readlength) {
  exit 1, "Read length not set, the provided value is '${params.readlength}'. Please specify a valid value for `--readlength`"
}

// Check if star_index is provided. (this is only when bam if not given)
if (!params.bams) {
  if (params.star_index) {
    star_index = params.star_index
  }else{
    exit 1, "STAR index path is required, Not provided. Please specify a valid value for `--star_index`"
  }
}else{ 
  star_index = false 
}

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

log.info "Splicing-pipelines - N F  ~  version 0.1"
log.info "====================================="
log.info "Run name                    : ${params.run_name}"
log.info "Date                        : ${date}"
log.info "Final prefix                : ${run_prefix}"
log.info "Assembly name               : ${params.assembly_name}"
log.info "Reads                       : ${params.reads}"
log.info "Bams                        : ${params.bams}"
log.info "Single-end                  : ${download_from('tcga') ? 'Will be checked for each TCGA BAM file' : params.singleEnd}"
log.info "GTF                         : ${params.gtf}"
log.info "STAR index                  : ${star_index}"
log.info "Stranded                    : ${params.stranded}"
log.info "rMATS pairs file            : ${params.rmats_pairs ? params.rmats_pairs : 'Not provided'}"
log.info "Adapter                     : ${download_from('tcga') ? 'Will be set for each sample based based on whether the sample is paired or single-end' : adapter_file}"
log.info "Read Length                 : ${params.readlength}"
log.info "Overhang                    : ${overhang}"
log.info "Minimum length              : ${minlen}"
log.info "Sliding window              : ${params.slidingwindow}"
log.info "rMATS variable read length  : ${variable_read_length}"
log.info "rMATS statoff               : ${params.statoff}"
log.info "rMATS paired stats          : ${params.paired_stats}"
log.info "rMATS novel splice sites    : ${params.novelSS}"
log.info "rMATS Minimum Intron Length : ${params.mil}"
log.info "rMATS Maximum Exon Length   : ${params.mel}"
log.info "Mismatch                    : ${params.mismatch}"
log.info "filterScore                 : ${params.filterScore}"
log.info "sjdbOverhangMin             : ${params.sjdbOverhangMin}"
log.info "STAR memory                 : ${params.star_memory ? star_memory : 'Not provided, Using STAR task max memory'}"
log.info "Test                        : ${params.test}"
log.info "Download from               : ${params.download_from ? params.download_from : 'FASTQs directly provided'}"
log.info "Key file                    : ${params.key_file ? params.key_file : 'Not provided'}"
log.info "Outdir                      : ${params.outdir}"
log.info "Max CPUs                    : ${params.max_cpus}"
log.info "Max memory                  : ${params.max_memory}"
log.info "Max time                    : ${params.max_time}"
log.info "Mega time                   : ${params.mega_time}"
log.info "Google Cloud disk-space     : ${params.gc_disk_size}"
log.info ""
log.info "\n"

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

if (params.download_from) {
  Channel
    .fromPath(params.reads)
    .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
    .splitCsv(skip:1)
    .map { sample -> sample[0].trim() }
    .set { accession_ids }
} 
// TODO: combine single and paired-end channel definitions
if (!params.download_from && params.singleEnd && !params.bams) {
  Channel
    .fromPath(params.reads)
    .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
    .splitCsv(skip:1)
    .map { sample_id, fastq -> [sample_id, file(fastq), params.singleEnd] }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
} 
if (!params.download_from && !params.singleEnd && !params.bams) {
  Channel
    .fromPath(params.reads)
    .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
    .splitCsv(skip:1)
    .map { sample_id, fastq1, fastq2 -> [ sample_id, [file(fastq1),file(fastq2)], params.singleEnd ] }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
}
if (params.bams) {
  Channel
    .fromPath(params.bams)
    .ifEmpty { exit 1, "Cannot find BAMs csv file : ${params.bams}" }
    .splitCsv(skip:1)
    .map { name, bam, bai -> [ name, file(bam), file(bai) ] }
    .into { indexed_bam; indexed_bam_rmats }
} 
Channel
  .from(params.assembly_name)
  .ifEmpty { exit 1, "Genome assembly name not set"}
  .set { assembly_name }
Channel
  .fromPath(params.gtf)
  .ifEmpty { exit 1, "Cannot find GTF file: ${params.gtf}" }
  .into { gtf_star ; gtf_stringtie; gtf_stringtie_merge; gtf_to_combine }
if (!params.bams) {
  Channel
    .fromPath(star_index)
    .ifEmpty { exit 1, "STAR index not found: ${star_index}" }
    .set { star_index }
}
Channel
  .fromPath(key_file)
  .ifEmpty { exit 1, "Key file not found: ${key_file}" }
  .set { key_file }
Channel
  .fromPath(params.multiqc_config)
  .ifEmpty { exit 1, "MultiQC config YAML file not found: ${params.multiqc_config}" }
  .set { multiqc_config }
if (params.rmats_pairs) {
  Channel
    .fromPath(params.rmats_pairs)
    .ifEmpty { exit 1, "Cannot find rMATS pairs file : ${params.rmats_pairs}" }
    .splitCsv(sep:' ')
    .map { row -> 
      def rmats_id = row[0]
      def b1 = row[1].toString().split(',')
      def b2 = row[2].toString().split(',')
      [ rmats_id, b1, b2 ]
    }
    .set { samples}
}

/*--------------------------------------------------
  Download FASTQs from GTEx or SRA
---------------------------------------------------*/

if ( download_from('gtex') || download_from('sra') ) {
  process get_accession {
    tag "${accession}"
    label 'tiny_memory'
    
    input:
    val(accession) from accession_ids
    each file(key_file) from key_file
    
    output:
    set val(accession), file(output_filename), val(params.singleEnd) into raw_reads_fastqc, raw_reads_trimmomatic

    script:
    def ngc_cmd_with_key_file = key_file.name != 'no_key_file.txt' ? "--ngc ${key_file}" : ''
    output_filename = params.singleEnd ? "${accession}.fastq.gz" : "${accession}_{1,2}.fastq.gz"
    """
    prefetch $ngc_cmd_with_key_file $accession --progress -o $accession
    fasterq-dump $ngc_cmd_with_key_file $accession --threads ${task.cpus} --split-3
    pigz *.fastq
    """
  }
}

/*--------------------------------------------------
  Download BAMs from TCGA
---------------------------------------------------*/

if (download_from('tcga')) {
  process get_tcga_bams {
    tag "${accession}"
    label 'low_memory'
    
    input:
    val(accession) from accession_ids
    each file(key_file) from key_file
    
    output:
    set val(accession), file("*.bam"), env(singleEnd) into bamtofastq
    file("${accession}_paired_info.csv") into paired_info

    script:
    // TODO: improve download speed by using `-n N_CONNECTIONS`
    // See https://github.com/IARCbioinfo/GDC-tricks#to-speed-up-the-download
    key_flag = key_file.name != 'no_key_file.txt' ? "-t $key_file" : ""
    """
    gdc-client download $accession $key_flag
    mv $accession/*.bam .

    # Check if reads are single or paired-end
    n_single_reads=\$(samtools view -c -F 1 *.bam)
    n_paired_reads=\$(samtools view -c -f 1 *.bam)

    singleEnd=true
    if (( \$n_paired_reads > \$n_single_reads )); then
        singleEnd=false
    fi

    echo "sample_id,n_single_reads,n_paired_reads,single_end" > ${accession}_paired_info.csv
    echo "$accession,\$n_single_reads,\$n_paired_reads,\$singleEnd" >> ${accession}_paired_info.csv
    """
  }

  paired_info
    .collectFile(name: "${params.outdir}/QC/tcga/paired_info.csv", keepHeader: true, skip: 1)
}

/*--------------------------------------------------
  Bedtools to extract FASTQ from BAM
---------------------------------------------------*/

if (download_from('tcga')) {
  process bamtofastq {
    tag "${name}"
    label 'low_memory'
    
    input:
    set val(name), file(bam), val(singleEnd) from bamtofastq
    
    output:
    set val(name), file("*.fastq.gz"), val(singleEnd) into raw_reads_fastqc, raw_reads_trimmomatic

    script:
    // 2GB reserved for Javaruntime and 2GB aside
    usable_mem = "${task.memory.toMega() - 4000}"
    // samtools takes memory per thread, so - 
    per_thread_mem = "${usable_mem.toInteger()/task.cpus.toInteger()}"
    // check end
    singleEnd=singleEnd.toBoolean()
    if (singleEnd) {
      """
      bedtools bamtofastq -i $bam -fq ${name}.fastq
      pigz *.fastq
      """
    } else {
      """
      samtools sort -@ ${task.cpus} -m ${per_thread_mem}M -n $bam > ${name}_sorted.bam
      bedtools bamtofastq \
        -i ${name}_sorted.bam \
        -fq ${name}_1.fastq \
        -fq2 ${name}_2.fastq
      pigz *.fastq
      """
    }
  }
}

if (!params.bams){

  /*--------------------------------------------------
    FastQC for quality control of input reads
  ---------------------------------------------------*/

  process fastqc {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/QC/raw", mode: 'copy'

    input:
    set val(name), file(reads), val(singleEnd) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results_raw

    script:
    """
    fastqc --casava --threads $task.cpus $reads
    """
  }

  /*--------------------------------------------------
    Trimmomatic to trim input reads
  ---------------------------------------------------*/

  raw_reads_trimmomatic
    .map { name, reads, singleEnd ->
      adapter = params.adapter ? file(params.adapter) : singleEnd ? file("$baseDir/adapters/TruSeq3-SE.fa") : file("$baseDir/adapters/TruSeq3-PE.fa")
      [ name, reads, singleEnd, adapter ]
    }
    .set { raw_reads_trimmomatic_adapter }

  process trimmomatic {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    set val(name), file(reads), val(singleEnd), file(adapter) from raw_reads_trimmomatic_adapter

    output:
    set val(name), file(output_filename), val(singleEnd) into (trimmed_reads_fastqc, trimmed_reads_star)
    file ("logs/${name}_trimmomatic.log") into trimmomatic_logs

    script:
    mode = singleEnd ? 'SE' : 'PE'
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
      MINLEN:${minlen} \
      CROP:${params.readlength}

    mkdir logs
    cp .command.log logs/${name}_trimmomatic.log
    """
  }

  /*--------------------------------------------------
    FastQC for quality control of input reads
  ---------------------------------------------------*/

  process fastqc_trimmed {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/QC/trimmed", mode: 'copy'

    input:
    set val(name), file(reads), val(singleEnd) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results_trimmed

    script:
    """
    fastqc --casava --threads $task.cpus $reads
    """
  }

  /*--------------------------------------------------
    STAR to align trimmed reads
  ---------------------------------------------------*/

  process star {
    tag "$name"
    label 'mega_memory'
    publishDir "${params.outdir}/star_mapped/${name}", mode: 'copy'
    publishDir "${params.outdir}/star_mapped/", mode: 'copy',
      saveAs: {filename -> 
          if (filename.indexOf(".bw") > 0) "all_bigwig/${name}.bw"
      }
    
    input:
    set val(name), file(reads), val(singleEnd) from trimmed_reads_star
    each file(index) from star_index
    each file(gtf) from gtf_star

    output:
    set val(name), file("${name}.Aligned.sortedByCoord.out.bam"), file("${name}.Aligned.sortedByCoord.out.bam.bai") into (indexed_bam, indexed_bam_rmats)
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    file "*Unmapped*" optional true
    file "${name}.bw"

    script:
    // TODO: check when to use `--outWigType wiggle` - for paired-end stranded stranded only?
    // TODO: find a better solution to needing to use `chmod`
    out_filter_intron_motifs = params.stranded ? '' : '--outFilterIntronMotifs RemoveNoncanonicalUnannotated'
    out_sam_strand_field = params.stranded ? '' : '--outSAMstrandField intronMotif'
    xs_tag_cmd = params.stranded ? "samtools view -h ${name}.Aligned.sortedByCoord.out.bam | gawk -v strType=2 -f /usr/local/bin/tagXSstrandedData.awk | samtools view -bS - > Aligned.XS.bam && mv Aligned.XS.bam ${name}.Aligned.sortedByCoord.out.bam" : ''
    // Set maximum available memory to be used by STAR to sort BAM files
    star_mem = params.star_memory ? params.star_memory : task.memory
    avail_mem_bam_sort = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 2000000000}" : ''
    """
    # Decompress STAR index if compressed
    if [[ $index == *.tar.gz ]]; then
      tar -xvzf $index
    fi

    STAR \
      --genomeDir ${index.toString().minus('.tar.gz')} \
      --readFilesIn $reads \
      --readMatesLengthsIn NotEqual \
      --outFileNamePrefix ${name}. \
      --runThreadN $task.cpus \
      --readFilesCommand zcat \
      --sjdbGTFfile $gtf \
      --sjdbOverhang $overhang \
      --alignSJDBoverhangMin $params.sjdbOverhangMin \
      --outFilterScoreMinOverLread $params.filterScore \
      --outFilterMatchNminOverLread $params.filterScore \
      --outFilterMismatchNmax $params.mismatch \
      --outFilterMultimapNmax 20 \
      --alignMatesGapMax 1000000 \
      --outSAMattributes All \
      --outSAMtype BAM SortedByCoordinate \
      $avail_mem_bam_sort \
      --outBAMsortingThreadN $task.cpus \
      --outFilterType BySJout \
      --twopassMode Basic \
      --alignEndsType EndToEnd \
      --alignIntronMax 1000000 \
      --outReadsUnmapped Fastx \
      --outWigType None $out_filter_intron_motifs $out_sam_strand_field

    chmod a+rw $name*
    $xs_tag_cmd
    samtools index ${name}.Aligned.sortedByCoord.out.bam
    bamCoverage -b ${name}.Aligned.sortedByCoord.out.bam -o ${name}.bw 
    """
  }
}

if (!params.test) {

  /*--------------------------------------------------
    Stringtie for transcript assembly and quantification 
  ---------------------------------------------------*/

  process stringtie {
    tag "$name"
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/${name}", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_index) from indexed_bam
    each file(gtf) from gtf_stringtie

    output:
    file "${name}.gtf" into stringtie_gtf
    file "${name}_for_DGE.gtf" into stringtie_dge_gtf

    script: 
    rf = params.stranded ? params.stranded == 'first-strand' ? '--rf' : '--fr' : ''
    """
    stringtie $bam -G $gtf -o ${name}.gtf $rf -a 8 -p $task.cpus
    stringtie $bam -G $gtf -o ${name}_for_DGE.gtf $rf -a 8 -e -p $task.cpus
    """
  }

  /*--------------------------------------------------
    Generate count matrix for all samples with prepDE.py
  ---------------------------------------------------*/

  process prep_de {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/count_matrix", mode: 'copy'

    input:
    file(gtf) from stringtie_dge_gtf.collect()

    output:
    file "sample_lst.txt"
    file "*gene_count_matrix.csv"
    file "*transcript_count_matrix.csv"

    script: 
    """
    echo "${gtf.join("\n").toString().replace("_for_DGE.gtf", "")}" > samples.txt
    echo "${gtf.join("\n")}" > gtfs.txt
    paste -d ' ' samples.txt gtfs.txt > sample_lst.txt
    prepDE.py -i sample_lst.txt  -l $params.readlength \
              -g ${run_prefix}_gene_count_matrix.csv -t ${run_prefix}_transcript_count_matrix.csv
    """
  } 

  /*--------------------------------------------------
    Stringtie merge GTF files
  ---------------------------------------------------*/

  process stringtie_merge {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/stringtie_merge", mode: 'copy'

    input:
    file('*.gtf') from stringtie_gtf.collect()
    file(gtf) from gtf_stringtie_merge

    output:
    file "gffcmp.annotated.corrected.gtf" into merged_gtf
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

  // Combine GTFs into a single channel so that rMATS runs twice (once for each GTF)
  gtf_rmats = gtf_to_combine.combine(merged_gtf).flatten()

  /*--------------------------------------------------
    rMATS to detect alternative splicing events
  ---------------------------------------------------*/

  if (params.rmats_pairs) {

    indexed_bam_rmats
      .map { name, bam, bai -> [name, bam] }
      .set { bam }
    
    // Group BAMs for each rMATS execution
    samples
      .map { row -> 
        def samples_rmats_id = []
        def rmats_id = row[0]
        def b1_samples = row[1]
        def b2_samples = row[2]
        b1_samples.each { sample ->
          samples_rmats_id.add([sample, 'b1', rmats_id])
        }
        b2_samples.each { sample ->
          samples_rmats_id.add([sample, 'b2', rmats_id])
        }
        samples_rmats_id
      }
      .flatMap()
      .combine(bam, by:0)
      .map { sample_id, b, rmats_id, bam -> [ rmats_id + b, rmats_id, bam] }
      .groupTuple()
      .map { b, rmats_id, bams -> [rmats_id[0], [b, bams]] }
      .groupTuple()
      .map { rmats_id, bams -> 
        def b1_bams = bams[0][0].toString().endsWith('b1') ? bams[0] : bams[1]
        def b2_bams = bams[0][0].toString().endsWith('b2') ? bams[0] : bams[1]
        rmats_id_bams = b2_bams == null ? [ rmats_id, b1_bams[1], "no b2", true ] : [ rmats_id, b1_bams[1] , b2_bams[1], false ]
        rmats_id_bams
      }
      .set { bams }

    process rmats {
      tag "$rmats_id ${gtf.simpleName}"
      label 'high_memory'
      publishDir "${params.outdir}/rMATS_out/${rmats_id}_${gtf.simpleName}", mode: 'copy'

      when:
      !params.skiprMATS

      input:
      set val(rmats_id), file(bams), file(b2_bams), val(b1_only) from bams
      each file(gtf) from gtf_rmats

      output:
      file "*.{txt,csv}" into rmats_out

      script:
      libType = params.stranded ? params.stranded == 'first-strand' ? 'fr-firststrand' : 'fr-secondstrand' : 'fr-unstranded'
      mode = params.singleEnd ? 'single' : 'paired'
      variable_read_length_flag = variable_read_length ? '--variable-read-length' : ''
      statoff = params.statoff ? '--statoff' : ''
      paired_stats = params.paired_stats ? '--paired-stats' : ''
      novelSS = params.novelSS ? '--novelSS' : ''   
      if (b1_only) {
        b1_bams = bams.join(",")
        b2_cmd = ''
        b2_flag = ''
        b2_config_cmd = ''
      } else {
        b1_bams = bams.join(",")
        b2_bams = b2_bams.join(",")
        b2_cmd = "echo $b2_bams > b2.txt"
        b2_flag = "--b2 b2.txt"
        b2_config_cmd = "echo b2 b2.txt >> \$rmats_config"
      }
      """
      echo $b1_bams > b1.txt
      $b2_cmd
      rmats.py \
        --b1 b1.txt $b2_flag \
        --gtf $gtf \
        --od ./ \
        --tmp tmp \
        --libType $libType \
        -t $mode \
        --nthread $task.cpus \
        --readLength ${params.readlength} \
        --mil ${params.mil} \
        --mel ${params.mel} $variable_read_length_flag $statoff $paired_stats $novelSS
      rmats_config="config_for_rmats_and_postprocessing.txt"
      echo b1 b1.txt > \$rmats_config
      $b2_config_cmd
      echo rmats_gtf       ${gtf} >> \$rmats_config
      echo ref_gtf         ${gtf} >> \$rmats_config
      echo fasta           ${params.assembly_name} >> \$rmats_config
      echo reads           ${params.singleEnd ? 'single' : 'paired'} >> \$rmats_config
      echo readlen         ${params.readlength} >> \$rmats_config
      echo rmats_id        ${rmats_id} >> \$rmats_config
      
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
      publishDir "${params.outdir}/rMATS_out/${name1}_vs_${name2}_${gtf.simpleName}", mode: 'copy'

      when:
      !params.skiprMATS

      input:
      set val(name1), file(bam1), val(name2), file(bam2) from paired_samples
      each file (gtf) from gtf_rmats

      output:
      file "*.{txt,csv}" into paired_rmats_out

      script:
      libType = params.stranded ? params.stranded == 'first-strand' ? 'fr-firststrand' : 'fr-secondstrand' : 'fr-unstranded'
      mode = params.singleEnd ? 'single' : 'paired'
      variable_read_length_flag = variable_read_length ? '--variable-read-length' : ''
      statoff = params.statoff ? '--statoff' : ''
      paired_stats = params.paired_stats ? '--paired-stats' : ''
      novelSS = params.novelSS ? '--novelSS' : ''   
      """
      ls $bam1 > b1.txt
      ls $bam2 > b2.txt
      rmats.py \
        --b1 b1.txt \
        --b2 b2.txt \
        --gtf $gtf \
        --od ./ \
        --tmp tmp \
        --libType $libType \
        -t $mode \
        --nthread $task.cpus \
        --readLength ${params.readlength} \
        --mil ${params.mil} \
        --mel ${params.mel} $variable_read_length_flag $statoff $paired_stats $novelSS
      rmats_config="config_for_rmats_and_postprocessing.txt"
      echo b1 b1.txt > \$rmats_config
      echo b2 b2.txt >> \$rmats_config
      echo rmats_gtf $gtf >> \$rmats_config
      echo rmats_gtf       ${gtf} >> \$rmats_config
      echo ref_gtf         ${gtf} >> \$rmats_config
      echo fasta           ${params.assembly_name} >> \$rmats_config
      echo reads           ${params.singleEnd ? 'single' : 'paired'} >> \$rmats_config
      echo readlen         ${params.readlength} >> \$rmats_config
      echo rmats_id        ${name1}_vs_${name2} >> \$rmats_config
      
      LU_postprocessing.R
      """
    }
  }
}

/*--------------------------------------------------
  MultiQC to generate a QC HTML report
---------------------------------------------------*/

if (!params.bams) {
  process multiqc {
    label 'mega_memory'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skipMultiQC

    input:
    file (fastqc:'fastqc/*') from fastqc_results_raw.collect().ifEmpty([])
    file (fastqc:'fastqc/*') from fastqc_results_trimmed.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file (multiqc_config) from multiqc_config

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc . --config $multiqc_config -m fastqc -m star
    cp multiqc_report.html ${run_prefix}_multiqc_report.html
    """
  }
}

// define helper function
def download_from(db) {
  download_from.toLowerCase().contains(db)
}
