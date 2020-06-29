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
      --run_name		    User specified name used as prefix for output files
				    (defaut: no prefix, only date and time)
      --download_from               Database to download FASTQ/BAMs from (available = 'TCGA', 'GTEX' or 'SRA', false) (string)
                                    (default: false)
      --key_file                    For downloading reads, use TCGA authentication token (TCGA) or dbGAP repository key (GTEx, path)
                                    (default: false)  
 
    Main arguments:
      --gtf                         Path to reference GTF file (path)
                                    (default: no gtf specified) 
      --assembly_name               Genome assembly name (available = 'GRCh38' or 'GRCm38', string)
                                    (default: false)
      --star_index                  Path to STAR index (path)
                                    (default: no index specified)
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
      --filterScore 		    Controls --outFilterScoreMinOverLread and outFilterMatchNminOverLread
				    (default: 0.66)
      --sjdOverhangMin		    Controls --alignSJDBoverhangMin (int)
				    (default: 8)

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

// Check if user has set adapter sequence. If not set is based on the value of the singleEnd parameter
adapter_file = params.adapter ? params.adapter : params.singleEnd ? "$baseDir/adapters/TruSeq3-SE.fa" : "$baseDir/adapters/TruSeq3-PE.fa"
overhang = params.overhang ? params.overhang : params.readlength - 1
download_from = params.download_from ? params.download_from : ""
key_file = params.key_file ? params.key_file : "$baseDir/examples/assets/no_key_file.txt"
minlen = params.minlen ? params.minlen : params.readlength
variable_read_length = minlen == params.readlength ? false : true

log.info "Splicing-pipelines - N F  ~  version 0.1"
log.info "====================================="
log.info "Assembly name               : ${params.assembly_name}"
log.info "Reads                       : ${params.reads}"
log.info "Bams                        : ${params.bams}"
log.info "Single-end                  : ${params.singleEnd}"
log.info "GTF                         : ${params.gtf}"
log.info "STAR index                  : ${params.star_index}"
log.info "Stranded                    : ${params.stranded}"
log.info "rMATS pairs file            : ${params.rmats_pairs ? params.rmats_pairs : 'Not provided'}"
log.info "Adapter                     : ${adapter_file}"
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
log.info "sjdOverhangMin              : ${params.sjdOverhangMin}"
log.info "Test                        : ${params.test}"
log.info "Download from               : ${params.download_from ? params.download_from : 'FASTQs directly provided'}"
log.info "Key file                    : ${params.key_file ? params.key_file : 'Not provided'}"
log.info "Outdir                      : ${params.outdir}"
log.info "Max CPUs                    : ${params.max_cpus}"
log.info "Max memory                  : ${params.max_memory}"
log.info "Max time                    : ${params.max_time}"
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
    .map { sample_id, fastq -> [sample_id, file(fastq)] }
    .into { raw_reads_fastqc; raw_reads_trimmomatic }
} 
if (!params.download_from && !params.singleEnd && !params.bams) {
  Channel
    .fromPath(params.reads)
    .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
    .splitCsv(skip:1)
    .map { sample_id, fastq1, fastq2 -> [ sample_id, [file(fastq1),file(fastq2)] ] }
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
  .fromPath(adapter_file)
  .ifEmpty { exit 1, "Cannot find adapter file: ${adapter_file}" }
  .set { adapter }
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
    .fromPath(params.star_index)
    .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
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
    
    input:
    val(accession) from accession_ids
    each file(key_file) from key_file
    
    output:
    set val(accession), file("*.fastq.gz") into raw_reads_fastqc, raw_reads_trimmomatic

    script:
    def vdbConfigCmd = key_file.name != 'no_key_file.txt' ? "vdb-config --import ${key_file} ./" : ''
    """
    $vdbConfigCmd
    fasterq-dump $accession --threads ${task.cpus} --split-3
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
    
    input:
    val(accession) from accession_ids
    each file(key_file) from key_file
    
    output:
    set val(accession), file("*.bam") into bamtofastq

    script:
    // TODO: improve download speed by using `-n N_CONNECTIONS`
    // See https://github.com/IARCbioinfo/GDC-tricks#to-speed-up-the-download
    key_flag = key_file.name != 'no_key_file.txt' ? "-t $key_file" : ""
    """
    gdc-client download $accession $key_flag
    mv $accession/*.bam .
    """
  }
}

/*--------------------------------------------------
  Bedtools to extract FASTQ from BAM
---------------------------------------------------*/

if (download_from('tcga')) {
  process bamtofastq {
    tag "${name}"
    
    input:
    set val(name), file(bam) from bamtofastq
    
    output:
    set val(name), file("*.fastq.gz") into raw_reads_fastqc, raw_reads_trimmomatic

    script:
    if (params.singleEnd) {
      """
      bedtools bamtofastq -i $bam -fq ${name}.fastq
      pigz *.fastq
      """
    } else {
      """
      samtools sort --threads ${task.cpus} -n $bam > ${name}_sorted.bam
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
    label 'process_medium'
    publishDir "${params.outdir}/QC/raw", mode: 'copy'

    input:
    set val(name), file(reads) from raw_reads_fastqc

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

  process trimmomatic {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    set val(name), file(reads) from raw_reads_trimmomatic
    each file(adapter) from adapter

    output:
    set val(name), file(output_filename) into (trimmed_reads_fastqc, trimmed_reads_star)
    file ("logs/${name}_trimmomatic.log") into trimmomatic_logs

    script:
    mode = params.singleEnd ? 'SE' : 'PE'
    out = params.singleEnd ? "${name}_trimmed.fastq.gz" : "${name}_trimmed_R1.fastq.gz ${name}_unpaired_R1.fastq.gz ${name}_trimmed_R2.fastq.gz ${name}_unpaired_R2.fastq.gz"
    output_filename = params.singleEnd ? "${name}_trimmed.fastq.gz" : "${name}_trimmed_R{1,2}.fastq.gz"
    slidingwindow = params.slidingwindow ? 'SLIDINGWINDOW:4:15' : ''
    """
    trimmomatic \
      $mode \
      -threads $task.cpus \
      -phred33 \
      $reads \
      $out \
      ILLUMINACLIP:${adapter}:2:30:10:8:true \
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
    label 'process_medium'
    publishDir "${params.outdir}/QC/trimmed", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_reads_fastqc

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
    label 'high_memory'
    publishDir "${params.outdir}/star_mapped/${name}", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_reads_star
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
    xs_tag_cmd = params.stranded ? "samtools view -h ${name}.Aligned.sortedByCoord.out.bam | awk -v strType=2 -f /usr/local/bin/tagXSstrandedData.awk | samtools view -bS - > Aligned.XS.bam && mv Aligned.XS.bam ${name}.Aligned.sortedByCoord.out.bam" : ''
    endsType = variable_read_length ? 'Local' : 'EndToEnd'
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
      --alignSJoverhangMin $params.sjdOverhangMin \
      --outFilterScoreMinOverLread $params.filterScore \
      --outFilterMatchNminOverLread $params.filterScore \
      --outFilterMismatchNmax $params.mismatch \
      --outFilterMultimapNmax 20 \
      --alignMatesGapMax 1000000 \
      --outSAMattributes All \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 100000000000 \
      --outBAMsortingThreadN $task.cpus \
      --outFilterType BySJout \
      --twopassMode Basic \
      --alignEndsType $endsType \
      --alignIntronMax 1000000 \
      --outReadsUnmapped Fastx \
      --outWigType wiggle $out_filter_intron_motifs $out_sam_strand_field

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
    label 'process_medium'
    publishDir "${params.outdir}/star_mapped/${name}", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_index) from indexed_bam
    each file(gtf) from gtf_stringtie

    output:
    file "${name}.gtf" into stringtie_gtf
    file "${name}_for_DGE.gtf" into stringtie_dge_gtf

    script: 
    rf = params.stranded ? params.stranded == 'first-strand' ? '--rf' : '--fr' : ''
    rf = params.stranded ? '--rf' : ''
    """
    stringtie $bam -G $gtf -o ${name}.gtf $rf -a 8 -p $task.cpus
    stringtie $bam -G $gtf -o ${name}_for_DGE.gtf $rf -a 8 -e -p $task.cpus
    """
  }

  /*--------------------------------------------------
    Generate count matrix for all samples with prepDE.py
  ---------------------------------------------------*/

  process prep_de {
    label 'process_medium'
    publishDir "${params.outdir}/star_mapped/count_matrix", mode: 'copy'

    input:
    file(gtf) from stringtie_dge_gtf.collect()

    output:
    file "sample_lst.txt"
    file "gene_count_matrix.csv"
    file "transcript_count_matrix.csv"

    script: 
    """
    echo "${gtf.join("\n").toString().replace("_for_DGE.gtf", "")}" > samples.txt
    echo "${gtf.join("\n")}" > gtfs.txt
    paste -d ' ' samples.txt gtfs.txt > sample_lst.txt
    prepDE.py -i sample_lst.txt  -l $params.readlength
    """
  } 

  /*--------------------------------------------------
    Stringtie merge GTF files
  ---------------------------------------------------*/

  process stringtie_merge {
    label 'process_medium'
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
        [ rmats_id, b1_bams[1] + b2_bams[1] ]
      }
      .set { bams }

    process rmats {
      label 'high_memory'
      publishDir "${params.outdir}/rMATS_out/${rmats_id}_${gtf.simpleName}", mode: 'copy'
      tag "$rmats_id ${gtf.simpleName}"

      when:
      !params.skiprMATS

      input:
      set val(rmats_id), file(bams) from bams
      each file(gtf) from gtf_rmats

      output:
      file "*"

      script:
      libType = params.stranded ? params.stranded == 'first-strand' ? 'fr-firststrand' : 'fr-secondstrand' : 'fr-unstranded'
      mode = params.singleEnd ? 'single' : 'paired'
      variable_read_length_flag = variable_read_length ? '--variable-read-length' : ''
      statoff = params.statoff ? '--statoff' : ''
      paired_stats = params.paired_stats ? '--paired-stats' : ''
      novelSS = params.novelSS ? '--novelSS' : ''   
      n_samples_replicates = bams.size()
      n_replicates = n_samples_replicates.intdiv(2)
      bam_groups = bams.collate(n_replicates)
      b1_bams = bam_groups[0].join(",")
      b2_bams = bam_groups[1].join(",")
      """
      echo $b1_bams > b1.txt
      echo $b2_bams > b2.txt
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
      file "*"

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
    """
  }
}

// define helper function
def download_from(db) {
  download_from.toLowerCase().contains(db)
}
