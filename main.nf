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
      --reads                       Path to reads.csv file, which specifies the sample_id and path to FASTQ files
                                    for each read or read pair (path).
                                    This file is used if starting at beginning of pipeline. It can be file paths,
                                    s3 links or ftp link.
                                    (default: no reads.csv)
      --bams                        Path to bams.csv file which specifies sample_id and path to BAM and BAM.bai
                                    files (path)
                                    If this file is provided, pipeline will start at Stringtie (and proceed through
                                    rMATS and post processing).
                                    (default: no bams.csv)
      --rmats_pairs                 Path to rmats_pairs.txt file containing b1 (and b2) samples names (path)
                                    (default: no rmats_pairs specified)
      --run_name                    User specified name used as prefix for output files
                                    (defaut: no prefix, only date and time)
      --download_from               Database to download FASTQ/BAMs from (available = 'TCGA', 'GTEX' or 'GEN3-DRS',
                                    'SRA', 'FTP') (string)
                                    false should be used to run local files on the HPC (Sumner).
                                    'TCGA' can also be used to download GDC data including HCMI data.
                                    (default: false)
      --key_file                    For downloading reads, use TCGA authentication token (TCGA) or dbGAP repository
                                    key (GTEx, path) or credentials.json file in case of 'GEN3-DRS'
                                    (default: false)

    Main arguments:
      --gtf                         Path to reference GTF file (path)
                                    (default: no gtf specified)
      --assembly_name               Genome assembly name (available = 'GRCh38' or 'GRCm38', string)
                                    (default: false)
      --star_index                  Path to STAR index (path)
                                    Star indices must be generated prior to run (with correct STAR version)
                                    (default: false)
      --singleEnd                   Specifies that the input is single-end reads (bool)
                                    This parameter also automatically establishes the path to the SE or PE adapters.
                                    For PE, set to false.
                                    (default: false)
      --stranded                    Specifies that the input is stranded ('first-strand', 'second-strand',
                                    false (aka unstranded))
                                    'first-strand' refers to RF/fr-firststrand in this pipeline.
                                    (default: 'first-strand')
      --readlength                  Read length - Note that all reads will be cropped to this length(int)
                                    (default: no read length specified)
      -profile                      Configuration profile to use. Can use multiple (comma separated, string)
                                    On sumner, this should be set in the main.pbs or as a command-line parameter.
                                    Profile can only be activated from the command line.
                                    Available: base, docker, sumner, test and more.

    Trimmomatic:
      --minlen                      Drop the read if it is below a specified length (int)
                                    Default parameters turn on --variable-readlength
                                    To crop all reads and turn off --variable-readlength, set minlen = readlength
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
                                    For TCGA values:
                                    https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
                                    (default: 0.66)
      --sjdbOverhangMin             Controls --alignSJDBoverhangMin (int)
                                    For TCGA values:
                                    https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
                                    (default: 3)
      --soft_clipping               Enables soft clipping (bool)
                                    If true, the STAR parameter will be --alignEndsType 'Local' and the rMATS parameter
                                    --allow-clipping will be added.
                                    If false, the STAR parameter will be --alignEndsType 'EndToEnd' and no rMATS
                                    parameter is added.
                                    NOTE: Soft Clipping will cause read lengths to be variable, so turn soft_clipping
                                    off if reads need to be same length. Variable read length parameter is turned on
                                    in rMATS when minlen does not equal readlength.
                                    (default: true)
      --save_unmapped               Save unmapped and partially mapped reads in separate file (bool)
                                    (default: false)
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
                                    To run the first half of the pipeline (through STAR), set test = true.
                                    (default: false)
      --max_cpus                    Maximum number of CPUs (int)
                                    (default: 72)
      --max_memory                  Maximum memory (memory unit)
                                    (default: 760.GB)
      --max_time                    Maximum time (time unit)
                                    (default: 72.h)
      --skiprMATS                   Skip rMATS (bool)
                                    (default: false)
      --skipMultiQC                 Skip MultiQC (bool)
                                    (default: false)
      --outdir                      The output directory where the results will be saved (string)
                                    On Sumner, this must be set in the main.pbs or via command line.
                                    NF_splicing_pipeline.config will not overwrite main.pbs.
                                    (default: <directory where you submit the job>/results)
      --mega_time                   Sets time limit for processes withLabel 'mega_memory' in the main.nf using the
                                    base.config (time unit)
                                    Make sure '#SBATCH -t' in 'main.pbs' is appropriately set if you are changing this parameter.
                                    (default: 20.h)
      --gc_disk_size                Only specific to google-cloud executor. Adds disk-space for few aggregative processes.
                                    (default: "200 GB" based on 100 samples. Simply add 2 x Number of Samples)
      --debug                       This option will enable echo of script execution into STDOUT with some additional
                                    resource information (such as machine type, memory, cpu and disk space)
                                    (default: false)
      --error_strategy              Mode of pipeline handling failed processes.
                                    Possible values: 'terminate', 'finish', 'ignore', 'retry'.
                                    Check nextflow documnetation for detailed descriptions of each mode:
                                    https://www.nextflow.io/docs/latest/process.html#process-page-error-strategy
                                    Set this parameter in the main.pbs, on the command line, or see NF_splicing_pipeline.config
                                    example (does not work like normal config param)
                                    This does not overwrited CloudOS config, which is set to:
                                    'errorStrategy = { task.exitStatus in [3,9,10,14,143,137,104,134,139] ? 'retry': 'ignore'}
                                    (default (non-cloudos): 'finish')
      --cleanup                     This option will enable nextflow work folder cleanup upon pipeline successfull
                                    completion. All intermediate files from nexftlow processes' workdirs will be
                                    cleared, staging folder with staged files will not be cleared.
                                    If pipeline is completed with errors or interrupted cleanup will not be executed.
                                    Following successfull run resumed from the failed run with --cleanup option enabled
                                    will only clear folders of processess created in the latest run, it will not clear
                                    cached folders coming from previous pipleine runs.
                                    Set this parameter in the main.pbs, on the command line, or see NF_splicing_pipeline.config
                                    example (does not work like normal config param)
                                    (default non-cloudos: true; cloudos: false)


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

// Check if error_strategy parameter has a correct value
if (!params.allowed_error_strategies.contains(params.error_strategy)) {
  exit 1, "Error strategy \"${params.error_strategy}\" is not correct. Please choose one of: ${params.allowed_error_strategies.join(", ")}."
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
if (params.stranded) {log.info "strType                     : ${params.strType[params.stranded].strType}"}
log.info "Soft_clipping               : ${params.soft_clipping}"
log.info "Save unmapped               : ${params.save_unmapped}"
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
log.info "Max Retries                 : ${params.max_retries}"
log.info "Max CPUs                    : ${params.max_cpus}"
log.info "Max memory                  : ${params.max_memory}"
log.info "Max time                    : ${params.max_time}"
log.info "Mega time                   : ${params.mega_time}"
log.info "Google Cloud disk-space     : ${params.gc_disk_size}"
log.info "Debug                       : ${params.debug}"
log.info "Error strategy              : ${config.process.errorStrategy}"
log.info "Workdir cleanup             : ${params.cleanup}"
log.info ""
log.info "\n"

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

if (params.download_from) {
  if(download_from('gtex') || download_from('sra') || download_from('tcga') ){
      Channel
        .fromPath(params.reads)
        .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
        .splitCsv(skip:1)
        .map { sample -> sample[0].trim() }
        .set { accession_ids }
  }
  if(download_from('gen3-drs')){
      Channel
        .fromPath(params.reads)
        .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
        .splitCsv(skip:1)
        .map { md5sum, file_name, obj_id, file_size -> [md5sum, file_name, obj_id, file_size] }
        .set { ch_gtex_gen3_ids }
  }
  if(download_from('ftp')){
    Channel
        .fromPath(params.reads)
        .ifEmpty { exit 1, "Cannot find CSV reads file : ${params.reads}" }
        .splitCsv(skip:1)
        .map { sample -> sample[0].trim() }
        .set { accession_ids_ftp }
  }
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

if ( download_from('gen3-drs')) {
    if(!params.genome_fasta){
    exit 1, "A genome fasta file must be provided in order to convert CRAM files in GEN3-DRS download step."
    }
    Channel
        .fromPath(params.genome_fasta)
        .ifEmpty { exit 1, "${params.genome_fasta} is not present" }
        .set {ch_genome_fasta}
}

if ( download_from('sra')) {
    Channel
        .value(file(params.sra_config_file))
        .set {ch_sra_config_file}
}


/*--------------------------------------------------
  Download FASTQs from GTEx or SRA
---------------------------------------------------*/

if ( download_from('gtex') || download_from('sra') ) {
  process get_accession {
    publishDir "${params.outdir}/process-logs/${task.process}/${accession}/", pattern: "command-logs-*", mode: 'copy'

    tag "${accession}"
    label 'tiny_memory'
    
    input:
    val(accession) from accession_ids
    each file(key_file) from key_file
    file(sra_config) from ch_sra_config_file
    
    output:
    set val(accession), file(output_filename), val(params.singleEnd) into raw_reads_fastqc, raw_reads_trimmomatic
	  file("command-logs-*") optional true

    script:
    def ngc_cmd_with_key_file = key_file.name != 'no_key_file.txt' ? "--ngc ${key_file}" : ''
    output_filename = params.singleEnd ? "${accession}.fastq.gz" : "${accession}_{1,2}.fastq.gz"
    """
    mkdir .ncbi
    mv ${sra_config} .ncbi/
    prefetch $ngc_cmd_with_key_file $accession --progress -o $accession
    fasterq-dump $ngc_cmd_with_key_file $accession --threads ${task.cpus} --split-3
    pigz *.fastq

    # save .command.* logs
    ${params.savescript}
    """
  }
} 

if ( download_from('ftp') ) {
  process get_ftp_accession {
    publishDir "${params.outdir}/process-logs/${task.process}/${accession}/", pattern: "command-logs-*", mode: 'copy'
    
    tag "${accession}"
    label 'tiny_memory'

    input:
    val(accession) from accession_ids_ftp

    output:
    set val(accession), file(output_filename), val(params.singleEnd) into raw_reads_fastqc, raw_reads_trimmomatic

    script:
    output_filename = params.singleEnd ? "${accession}.fastq.gz" : "${accession}_{1,2}.fastq.gz"
    isSingle = params.singleEnd ? "true" : "false"

    """
    PREFIX="\$(echo "$accession" | head -c 6)"
    FTP_PATH="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/\${PREFIX}"
    SAMPLE=$accession

    if [[ "\${#SAMPLE}" == "9" ]]; then
      FTP_PATH="\${FTP_PATH}/$accession/$accession"
    elif [[ "\${#SAMPLE}" == "10" ]]; then
      SUFFIX="\${SAMPLE: -1}" 
      FTP_PATH="\${FTP_PATH}/00\${SUFFIX}/$accession/$accession"
    elif [[ "\${#SAMPLE}" == "11" ]]; then
      SUFFIX="\${SAMPLE: -2}" 
      FTP_PATH="\${FTP_PATH}/0\${SUFFIX}/$accession/$accession"      
    else    
      SUFFIX="\${SAMPLE: -3}"    
      FTP_PATH="\${FTP_PATH}/\${SUFFIX}/$accession/$accession"
    fi

    echo \$SAMPLE
    echo \$FTP_PATH

    if [ "$isSingle" = true ] ; then
      {
        wget "\${FTP_PATH}.fastq.gz"
      } || {
        wget "\${FTP_PATH}_1.fastq.gz" 
        mv ${accession}_1.fastq.gz ${accession}.fastq.gz
      }
    else
      wget "\${FTP_PATH}_1.fastq.gz"
      wget "\${FTP_PATH}_2.fastq.gz"
    fi
    """
  }
}

/*--------------------------------------------------
  Download BAMs from GTEx using GEN3_DRS 
---------------------------------------------------*/

if ( download_from('gen3-drs')) {
  process gen3_drs_fasp {
      tag "${file_name}"
      label 'low_memory'
      publishDir "${params.outdir}/process-logs/${task.process}/${file(file_name).baseName}", pattern: "command-logs-*", mode: 'copy'

      input:
      set val(md5sum), val(file_name), val(obj_id), val(obj_id), val(file_size) from ch_gtex_gen3_ids
      each file(key_file) from key_file
      each file(genome_fasta) from ch_genome_fasta
      
      output:
      set stdout, file("*.bam"), val(false) into bamtofastq
      file("command-logs-*") optional true
      
      script:
      """
      sample_name=\$(echo ${file_name} | cut -f1 -d".")
      
      drs_url=\$(python /fasp-scripts/fasp/scripts/get_drs_url.py ${obj_id} gcp_id ${key_file})
      signed_url=\$(echo \$drs_url | awk '\$1="";1')
      
      if [[ \$signed_url == *".bam"* ]]; then
          wget -O \${sample_name}.bam \$(echo \$signed_url)
          file_md5sum=\$(md5sum \${sample_name}.bam)
          if [[ ! "\$file_md5sum" =~ ${md5sum} ]]; then echo "md5 sum verification failed" > md5sum_check.log; exit 1; else echo "file is good" > md5sum_check.log; fi
      fi
      
      if [[ \$signed_url == *".cram"* ]]; then
          wget -O \${sample_name}.cram \$(echo \$signed_url)
          file_md5sum=\$(md5sum \${sample_name}.cram)
          if [[ ! "\$file_md5sum" =~ ${md5sum} ]]; then exit 1; else echo "file is good"; fi
          samtools view -b -T ${genome_fasta} -o \${sample_name}.bam \${sample_name}.cram
      fi

      # save .command.* logs
      ${params.savescript}

      printf "\$sample_name"
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
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'
    
    input:
    val(accession) from accession_ids
    each file(key_file) from key_file
    
    output:
    set val(accession), file("*.bam"), env(singleEnd) into bamtofastq
    file("${accession}_paired_info.csv") into paired_info
	  file("command-logs-*") optional true

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

    # save .command.* logs
    ${params.savescript}
    """
  }

  paired_info
    .collectFile(name: "${params.outdir}/QC/tcga/paired_info.csv", keepHeader: true, skip: 1)
}

/*--------------------------------------------------
  Bedtools to extract FASTQ from BAM
---------------------------------------------------*/

if (download_from('tcga') || download_from('gen3-drs')) {
  process bamtofastq {
    tag "${name}"
    label 'mid_memory'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}/", pattern: "command-logs-*", mode: 'copy'

    input:
    set val(name), file(bam), val(singleEnd) from bamtofastq
    
    output:
    set val(name), file("*.fastq.gz"), val(singleEnd) into raw_reads_fastqc, raw_reads_trimmomatic
	  file("command-logs-*") optional true

    script:
    // samtools takes memory per thread
    // 6GB reserved for Javaruntime
    def usable_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
    def per_thread_mem = (task.memory && usable_mem) ? "${(task.memory.toBytes() - 6000000000) / (task.cpus * 2)}" : ''
    // check end
    singleEnd=singleEnd.toBoolean()
    if (singleEnd) {
      """
      bedtools bamtofastq -i $bam -fq ${name}.fastq
      pigz *.fastq

      # save .command.* logs
      ${params.savescript}
      """
    } else {
      """
      samtools sort -@ ${task.cpus} -m ${per_thread_mem} -n $bam > ${name}_sorted.bam
      bedtools bamtofastq \
        -i ${name}_sorted.bam \
        -fq ${name}_1.fastq \
        -fq2 ${name}_2.fastq
      pigz *.fastq

      # save .command.* logs
      ${params.savescript}
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
    publishDir "${params.outdir}/QC/raw", pattern: "*_fastqc.{zip,html}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    input:
    set val(name), file(reads), val(singleEnd) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results_raw
    file("command-logs-*") optional true

    script:
    """
    fastqc --casava --threads $task.cpus $reads

    # save .command.* logs
    ${params.savescript}
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
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    input:
    set val(name), file(reads), val(singleEnd), file(adapter) from raw_reads_trimmomatic_adapter

    output:
    set val(name), file(output_filename), val(singleEnd) into (trimmed_reads_fastqc, trimmed_reads_star)
    file ("logs/${name}_trimmomatic.log") into trimmomatic_logs
    file("command-logs-*") optional true

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

    # save .command.* logs
    ${params.savescript}
    """
  }

  /*--------------------------------------------------
    FastQC for quality control of input reads
  ---------------------------------------------------*/

  process fastqc_trimmed {
    tag "$name"
    label 'low_memory'
    publishDir "${params.outdir}/QC/trimmed", pattern: "*_fastqc.{zip,html}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    input:
    set val(name), file(reads), val(singleEnd) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results_trimmed
    file("command-logs-*") optional true

    script:
    """
    fastqc --casava --threads $task.cpus $reads

    # save .command.* logs
    ${params.savescript}
    """
  }

  /*--------------------------------------------------
    STAR to align trimmed reads
  ---------------------------------------------------*/

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
    publishDir "${params.outdir}/star_mapped/${name}", pattern: "*{out.bam,out.bam.bai,out,ReadsPerGene.out.tab,SJ.out.tab}*" , mode: 'copy'
    publishDir "${params.outdir}/star_mapped/${name}", pattern: "*Unmapped*", mode: 'copy'
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
    file "*ReadsPerGene.out.tab"
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    file "*Unmapped*" optional true
    file "${name}.bw"
    file("command-logs-*") optional true

    script:
    // TODO: check when to use `--outWigType wiggle` - for paired-end stranded stranded only?
    // TODO: find a better solution to needing to use `chmod`
    out_filter_intron_motifs = params.stranded ? '' : '--outFilterIntronMotifs RemoveNoncanonicalUnannotated'
    out_sam_strand_field = params.stranded ? '' : '--outSAMstrandField intronMotif'
    xs_tag_cmd = params.stranded ? "samtools view -h ${name}.Aligned.sortedByCoord.out.bam | gawk -v q=${params.strType[params.stranded].strType} -f /usr/local/bin/tagXSstrandedData.awk | samtools view -bS - > Aligned.XS.bam && mv Aligned.XS.bam ${name}.Aligned.sortedByCoord.out.bam" : ''
    endsType = params.soft_clipping ? 'Local' : 'EndToEnd'
    // Set maximum available memory to be used by STAR to sort BAM files
    star_mem = params.star_memory ? params.star_memory : task.memory
    avail_mem_bam_sort = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 2000000000}" : ''
    save_unmapped_reads = params.save_unmapped ? '--outReadsUnmapped Fastx' : ''
    """
    ${pre_script_run_resource_status}

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
      --alignEndsType $endsType \
      --alignIntronMax 1000000 \
      $save_unmapped_reads \
      --quantMode GeneCounts \
      --outWigType None $out_filter_intron_motifs $out_sam_strand_field

    chmod a+rw $name*
    $xs_tag_cmd
    samtools index ${name}.Aligned.sortedByCoord.out.bam
    bamCoverage -b ${name}.Aligned.sortedByCoord.out.bam -o ${name}.bw 

    ${post_script_run_resource_status}
    rm -r ${file(index).name.minus('.gz').minus('.tar')}   # not simpleName or twice baseName because index has dot's in name:  star_2.7.9a_yeast_chr_I.tar.gz

    # save .command.* logs
    ${params.savescript}
    """
  }
}

if (!params.test) {

  /*--------------------------------------------------
    Stringtie for transcript assembly and quantification 
  ---------------------------------------------------*/

  process stringtie {
    tag "$name"
    label 'mega_memory'
    publishDir "${params.outdir}/star_mapped/${name}", pattern: "[!command-logs-]*", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/${name}", pattern: "command-logs-*", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_index) from indexed_bam
    each file(gtf) from gtf_stringtie

    output:
    file "${name}.gtf" into stringtie_gtf
    file "${name}_for_DGE.gtf" into stringtie_dge_gtf
    file("command-logs-*") optional true

    script: 
    rf = params.stranded ? params.stranded == 'first-strand' ? '--rf' : '--fr' : ''
    """
    stringtie $bam -G $gtf -o ${name}.gtf $rf -a 8 -p $task.cpus
    stringtie $bam -G $gtf -o ${name}_for_DGE.gtf $rf -a 8 -e -p $task.cpus

    # save .command.* logs
    ${params.savescript}
    """
  }

  /*--------------------------------------------------
    Generate count matrix for all samples with prepDE.py
  ---------------------------------------------------*/

  process prep_de {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/count_matrix", pattern: "{sample_lst.txt,*gene_count_matrix.csv,*transcript_count_matrix.csv}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    input:
    file(gtf) from stringtie_dge_gtf.collect()

    output:
    file "sample_lst.txt"
    file "*gene_count_matrix.csv"
    file "*transcript_count_matrix.csv"
    file("command-logs-*") optional true

    script: 
    """
    echo "${gtf.join("\n").toString().replace("_for_DGE.gtf", "")}" > samples.txt
    echo "${gtf.join("\n")}" > gtfs.txt
    paste -d ' ' samples.txt gtfs.txt > sample_lst.txt
    prepDE.py -i sample_lst.txt  -l $params.readlength \
              -g ${run_prefix}_gene_count_matrix.csv -t ${run_prefix}_transcript_count_matrix.csv

    # save .command.* logs
    ${params.savescript}
    """
  } 

  /*--------------------------------------------------
    Stringtie merge GTF files
  ---------------------------------------------------*/

  process stringtie_merge {
    label 'mid_memory'
    publishDir "${params.outdir}/star_mapped/stringtie_merge", pattern: "{gffcmp.annotated.corrected.gtf,gffcmp.*}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    input:
    file('*.gtf') from stringtie_gtf.collect()
    file(gtf) from gtf_stringtie_merge

    output:
    file "gffcmp.annotated.corrected.gtf" into merged_gtf
    file "gffcmp.*" into gffcmp
    file("command-logs-*") optional true

    script:
    """
    ls -1 *.gtf > assembly_gtf_list.txt
    stringtie --merge -G $gtf -o stringtie_merged.gtf assembly_gtf_list.txt -p $task.cpus
    gffcompare -R -V -r $gtf stringtie_merged.gtf
    correct_gene_names.R
    gffread -E gffcmp.annotated.corrected.gff -T -o gffcmp.annotated.corrected.gtf

    # save .command.* logs
    ${params.savescript}
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
      publishDir "${params.outdir}/rMATS_out/${rmats_id}_${gtf.simpleName}", pattern: "{*.txt,*.csv,tmp/*_read_outcomes_by_bam.txt}", mode: 'copy'
      publishDir "${params.outdir}/process-logs/${task.process}/${rmats_id}_${gtf.simpleName}", pattern: "command-logs-*", mode: 'copy'

      when:
      !params.skiprMATS

      input:
      set val(rmats_id), file(bams), file(b2_bams), val(b1_only) from bams
      each file(gtf) from gtf_rmats

      output:
      file "*.{txt,csv}" into rmats_out
      file "tmp/*_read_outcomes_by_bam.txt"
      file("command-logs-*") optional true

      script:
      libType = params.stranded ? params.stranded == 'first-strand' ? 'fr-firststrand' : 'fr-secondstrand' : 'fr-unstranded'
      mode = params.singleEnd ? 'single' : 'paired'
      variable_read_length_flag = variable_read_length ? '--variable-read-length' : ''
      statoff = params.statoff ? '--statoff' : ''
      paired_stats = params.paired_stats ? '--paired-stats' : ''
      novelSS = params.novelSS ? '--novelSS' : ''   
      allow_clipping = params.soft_clipping ? '--allow-clipping' : ''
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
        --mel ${params.mel} $variable_read_length_flag $statoff $paired_stats $novelSS $allow_clipping
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

      # save .command.* logs
      ${params.savescript}
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
      publishDir "${params.outdir}/rMATS_out/${name1}_vs_${name2}_${gtf.simpleName}", pattern: "{*.txt,*.csv,tmp/*_read_outcomes_by_bam.txt}", mode: 'copy'
      publishDir "${params.outdir}/process-logs/${task.process}/${name1}_vs_${name2}_${gtf.simpleName}", pattern: "command-logs-*", mode: 'copy'

      when:
      !params.skiprMATS

      input:
      set val(name1), file(bam1), val(name2), file(bam2) from paired_samples
      each file (gtf) from gtf_rmats

      output:
      file "*.{txt,csv}" into paired_rmats_out
      file "tmp/*_read_outcomes_by_bam.txt"
      file("command-logs-*") optional true

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

      # save .command.* logs
      ${params.savescript}
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
    publishDir "${params.outdir}/MultiQC", pattern: "{*multiqc_report.html,*_data/*,trimmomatic}", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    when:
    !params.skipMultiQC

    input:
    file (fastqc:'fastqc/*') from fastqc_results_raw.collect().ifEmpty([])
    file (fastqc:'fastqc/*') from fastqc_results_trimmed.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file (multiqc_config) from multiqc_config
    file ('trimmomatic/*') from trimmomatic_logs.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data/*"
    file ('trimmomatic')
    file("command-logs-*") optional true

    script:
    """
    multiqc . --config $multiqc_config -m fastqc -m star -m trimmomatic
    cp multiqc_report.html ${run_prefix}_multiqc_report.html

    # save .command.* logs
    ${params.savescript}
    """
  }
}


// Get tool versions

process collect_tool_versions_env1 {
    // TODO: This collects tool versions for only one base enviroment/container - 'gcr.io/nextflow-250616/splicing-pipelines-nf:gawk'
    // need to get tool versions from other enviroment/container
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    output:
    file("tool_versions.txt") into ch_tool_versions
    file("command-logs-*") optional true

    script:
    """
    touch tool_versions.txt
    conda list -n splicing-pipelines-nf | grep fastqc | tail -n 1 >> tool_versions.txt
    conda list -n splicing-pipelines-nf | grep trimmomatic | tail -n 1 >> tool_versions.txt
    conda list -n splicing-pipelines-nf | grep star | tail -n 1 >> tool_versions.txt
    conda list -n splicing-pipelines-nf | grep samtools | tail -n 1 >> tool_versions.txt
    conda list -n splicing-pipelines-nf | grep deeptools | tail -n 1 >> tool_versions.txt
    conda list -n splicing-pipelines-nf | grep multiqc | tail -n 1 >> tool_versions.txt
    conda list -n splicing-pipelines-nf | grep gffread | tail -n 1 >> tool_versions.txt
    echo -e "stringtie" ' \t\t\t\t ' \$(stringtie --version) >> tool_versions.txt

    # save .command.* logs
    ${params.savescript}
    """
}

process collect_tool_versions_env2 {
    echo true
    publishDir "${params.outdir}/tool-versions/env2/", pattern: "[!command-logs-]*", mode: 'copy'
    publishDir "${params.outdir}/process-logs/${task.process}/", pattern: "command-logs-*", mode: 'copy'

    input:
    file(tool_versions) from ch_tool_versions

    output:
    file("tool_versions.txt") into ch_all_tool_versions
    file("command-logs-*") optional true

    script:
    """
    conda list -n rmats4 | grep rmats | tail -n 1 >> tool_versions.txt

    # save .command.* logs
    ${params.savescript}
    """
}

// define helper function
def download_from(db) {
  download_from.toLowerCase().contains(db)
}


// Completion notification

workflow.onComplete {

    c_green = "\033[0;32m";
    c_purple = "\033[0;35m";
    c_red = "\033[0;31m";
    c_reset = "\033[0m";

    if (workflow.success) {
        log.info "-${c_purple}[splicing-pipelines-nf]${c_green} Pipeline completed successfully${c_reset}-"
        if (params.cleanup) {
          log.info "-${c_purple}[splicing-pipelines-nf]${c_green} Cleanup: Working directory cleared from intermediate files generated with current run: '${workflow.workDir}'  ${c_reset}-"
        }
    } else { // To be shown requires errorStrategy = 'finish'
        log.info "-${c_purple}[splicing-pipelines-nf]${c_red} Pipeline completed with errors${c_reset}-"
        if (params.cleanup) {
          log.info "-${c_purple}[splicing-pipelines-nf]${c_red} Cleanup: Working directory was not cleared from intermediate files due to pipeline errors. You can re-use them with -resume option.  ${c_reset}-"
        }
    }
}
