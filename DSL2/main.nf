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

    include { fastqc } from './modules/fastqc' //FastQC for quality control of input reads
    include { trimmomatic } from './modules/trimmomatic'
    include { fastqc_trimmed } from './modules/fastqc_trimmed'
    include { star } from './modules/star'
    include { index } from './modules/index'
    include { bamCoverage } from './modules/bamCoverage'
    include { stringtie } from './modules/stringtie'
    include { prep_de } from './modules/prep_de'
    include { stringtie_merge } from './modules/stringtie_merge'
    include { gffcompare } from './modules/gffcompare'
    include { R } from './modules/R'
    include { gffread } from './modules/gffread'

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
