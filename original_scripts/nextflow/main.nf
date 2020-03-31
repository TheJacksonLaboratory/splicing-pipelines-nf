/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'splicing-pipelines-nf' a pipeline repository to run Olga Anczukow's splicing pipeline.
 *
 * Main rmats-NF pipeline script
 *
 * @authors
 * Marina Yurieva <marina.yurieva@jax.org>
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 */

log.info "Splicing-pipelines - N F  ~  version 0.1"
log.info "====================================="
log.info "Accession list        : ${params.accessionList}"
log.info "Key file              : ${params.keyFile ? params.keyFile : 'Not provided'}"
log.info "Genome                : ${params.genomeFile}"
log.info "Read type             : ${params.readType}"
log.info "Read length           : ${params.readLength}"
log.info "output                : ${params.output}"
log.info "\n"

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run jacksonlabs/splicing-pipelines-nf --accessionList 'accession_list.txt' -profile docker
    Mandatory arguments:
      --accessionList               Path to input file with accession list to fetch from SRA
      --genomeFile                  Path to input genome file (version and url?)
      --keyFile                     Path to a keyfile used to fetch restricted access datasets with SRAtools
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: short-test, key-test, ...
    Generic:
      --skipTrimming                Sets if trimming has to be skipped
      --readType                    Specifies what type of input reads are to be processed: single end or paired end
      --readLength                  Specifies the read length
    """.stripIndent()
}

/*********************************
 *      CHANNELS SETUP           *
 *********************************/

alternativeSplicingTypeList = ['a3ss', 'a5ss', 'mxe', 'ri', 'se']
junctionCountTypeList       = ['jc','jcec']
countingTypeList            = ['ijc','sjc','inc','inclen','skiplen']

key_file = file(params.keyFile)

if(!params.starIndexPath) {
  exit 1, "Cannot find STAR index path for rMATs"
}

if(!params.splitNumber) {
    exit 1, "Cannot find a splitNumber value set for createMatrices"
}

if(params.accessionList) {
    Channel
        .fromPath( params.accessionList )
        .splitText()
        .map{ it.trim() }
        .dump(tag: 'AccessionList content')
        .set { accessionIDs }
} else {
    exit 1, "Accession list file not provided"
}

if(params.gencodeFile) {
    Channel
        .value(file(params.gencodeFile))
        .ifEmpty { error "Cannot find any Gencode GTF annotaton for parameter --gencode: ${params.gencodeFile}" }
        .set { gencodeFile }    
} else {
    exit 1, "Gencode annotation file not provided"
}

if (!params.readLength) {
    exit 1, "Read length parameter not provided"
}
if (!params.readType) {
    exit 1, "Read type parameter not provided"
}

/**********************************
 *      PIPELINE PROCESSES        *
 **********************************/

/*
 * Get accession samples from SRA with or without keyFile
 */
if(!params.reads) {
    process getAccession {
        tag "${accession}"
        
        input:
        val accession from accessionIDs
        file keyFile from key_file
        
        output:
        set val(accession), file("*.fastq.gz") into readFiles
        
        script:
        def vdbConfigCmd = keyFile.name != 'NO_FILE' ? "vdb-config --import ${keyFile} ./" : ''
        """
        $vdbConfigCmd
        fasterq-dump $accession --threads ${task.cpus} --split-3
        pigz *.fastq
        """
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.readType == 'single' ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set { readFiles }
}

/*
 * Optionally trim files to a fixed length with Fastp
 */
if(!params.skipTrimming) {
    process trimming {
        tag "${accession}"
        
        input:
        set val(accession), file(reads) from readFiles
        
        output:
        set val(accession), file("*_trimmed.fastq.gz") into trimmedFiles
        
        script:
        if (params.readType == "single") {
          """
          fastp -w ${task.cpus} -i ${reads} -b ${params.readLength} -o ${accession}_trimmed.fastq.gz
          """
        } else {
          """
          fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} -b ${params.readLength} -B  ${params.readLength} -o ${accession}_R1_trimmed.fastq.gz -O ${accession}_R2_trimmed.fastq.gz
          """
        }
    }
} else {
   readFiles
       .set {trimmedFiles}
}

/*
 * Map read files to the transcriptome with Hisat2
 */

process mapping {
    tag "mapping: $reads"
    
    input:
    set val(name), file(reads) from trimmedFiles
    
    output:
    set val(name), file("${name}.bam") into hisat2Bams
    file("${name}.hisat2_summary.txt") into hisat2Multiqc
    
    script:
    if (params.readType == "single") {
                """
            hisat2 -x ${params.hisat2index} \
                   -U ${reads} \
                   -p ${task.cpus} \
                   --met-stderr \
                   --new-summary \
                   --summary-file ${name}.hisat2_summary.txt \
                   | samtools view -bS - > ${name}.bam
                """
    } else {
        """
            hisat2 -x ${params.hisat2index} \
                   -1 ${reads[0]} \
                   -2 ${reads[1]} \
                   -p ${task.cpus} \
                   --met-stderr \
                   --new-summary \
                   --summary-file ${name}.hisat2_summary.txt \
                   | samtools view -bS - > ${name}.bam
        """
    }
}


/*
 * Sort BAM files with SAMTOOLS
 */

process sortbam {
    tag "sortbam: $name"
    
    input:
    set val(name), file(bam) from hisat2Bams
    
    output:
    set val(name), file("${name}.sorted.bam") into sortedBams
    
    script:
    avail_mem=""
    """
    samtools sort \\
        $bam \\
        -@ ${task.cpus} ${avail_mem} \\
        -o ${name}.sorted.bam
    samtools index ${name}.sorted.bam
    """
}

/*
 * Mark duplicates and sort with picard tools
 */

process markduplicates {
    tag "markdups: $name"
    
    input:
    set val(name), file(bam) from sortedBams
    
    output:
    set val(name), file("${name}.sorted.nodup.bam") into markedDups, bamstatsSamples
    file("${name}.metric.txt") into markduplicatesMultiqc
    
    script:
    // Runtime java Mark duplicates options
    markdup_java_options="-Xmx${task.memory.toGiga()}g"
    """
    picard ${markdup_java_options} MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${name}.sorted.nodup.bam \\
        METRICS_FILE=${name}.metric.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

/*
 * Generate BAM file statistics
 */
 
process bamstats {
    tag "bamstats: $name"

    publishDir = [path: "${params.output}/flagstats", mode: 'copy', overwrite: 'true' ]
    
    input:
    set val(name), file(bam) from bamstatsSamples
    
    output:
    file("${name}.sorted.nodup.bam.flagstat.txt") into flagstatMultiqc
    
    script:
    """
    samtools flagstat $bam > ${name}.sorted.nodup.bam.flagstat.txt
    samtools view $bam | cut -f 10 | perl -ne 'chomp;print length(\$_) . "\n"' | sort | uniq -c >> ${name}.sorted.nodup.bam.flagstat.txt
    """
}

/*
 * Group samples and their BAM files in pairs
 */

markedDups
     .toSortedList { entry -> entry[0] }
     .flatten()
     .collate(4, false)
     .set { pairedSamples }

/*
 * Run rMATS in pairs of samples
 */

process paired_rmats {
    tag "paired_rmats: ${sample1Name}_${sample2Name}"
    label 'rmats'

    publishDir = [path: "${params.output}/paired_rmats", mode: 'copy']

    input:
    set val(sample1Name), file(sample1Bam), val(sample2Name), file(sample2Bam) from pairedSamples
    file(gencodeGtf) from gencodeFile

    output:
    set val(sample1Name), val(sample2Name), file('rmats_output') into rmatsCounts
    file 'rmats_output/ASEvents/fromGTF*' into fromGtf

    script:
    """
    RNASeq-MATS.py \
        -b1 $sample1Bam \
        -b2 $sample2Bam \
        -gtf $gencodeGtf \
        -t ${params.readType} \
        -len ${params.readLength} \
        -c ${params.cutoffSplicingDifference} \
        -analysis ${params.analysisType} \
        -novelSS ${params.novelDetectionFlag} \
        -a ${params.anchorLength} \
        -keepTemp \
        -o rmats_output
    """
}
/*
 * Save counts per sample
 */

process sampleCountsSave {
    tag "sampleCountsSave: ${sample1Name}_${sample2Name}"
    label 'postrmats'

    publishDir = [path: "${params.output}/split_matrices", mode: 'copy', overwrite: 'true' ]

    input:
    set val(sample1Name), val(sample2Name), file(counts) from rmatsCounts
    
    output:
    set file("rmats_output/${sample1Name}.*.txt"), file("rmats_output/${sample2Name}.*.txt") into savedSampleCounts
    
    script:
    """
    sampleCountsSave.sh rmats_output ${sample1Name} ${sample2Name}
    """
}


/* 
 * Normalize the counts
 *          &
 * Create matrices from all sample count files
 */

 process createMatrices {
    tag "createMatrices: ${alternativeSplicingType}/${junctionCountType}/${countingType}/${params.splitNumber}"
    label 'postrmats'

    publishDir = [path: "${params.output}/merged_matrices", mode: 'copy', overwrite: 'true' ]

    input:
    file(allSamplesCounts) from savedSampleCounts.flatten().collect()
    file(fromGtfs) from fromGtf.first()
    each alternativeSplicingType from alternativeSplicingTypeList
    each junctionCountType from junctionCountTypeList
    each countingType from countingTypeList
    
    output:
    file 'rmats_final*' into splicingMatrices
    
    script:
    """
    create_matrices_from_files.sh ${alternativeSplicingType} ${junctionCountType} ${countingType} ${params.splitNumber}
    """
}
