#!/usr/bin/env nextflow
//rMATS_pipeline_samtools.nf
path2trim = file(params.path2trim)
params.reads = "PE"
adapters = file(params.adapters)
genomeDir = file(params.genomeDir)
genomeDirchr = file(params.genomeDirchr)
refgtf = file(params.refgtf)
readlength = "100"
params.bam2fastqlist = ''
params.fastqlist = 'testfastqlist.csv'
params.outgtf = 'SRproteins_MCF10A_novel.gtf'
outgtf = params.outgtf
if(params.bam2fastqlist) {
        Channel.fromPath(params.bam2fastqlist)
        .splitCsv()
        .map { line ->
                [file(line[0]), file(line[1])]
        }.set { bam_entries }
        process bam2fastq {
                container 'gcr.io/jax-awilliams-ctrl-32/star'
                cpus = 4
                publishDir "fastq", pattern:'*.*', mode:'copy'
                input:
                set file(bamfile), val(output_name) from bam_entries
                output:
                file "*.fastq.gz" into fastq
                script:
                """
                samtools fastq -@ ${task.cpus} -1 ${output_name}_R1.fastq.gz -2 ${output_name}_R2.fastq.gz $bamfile
                """
                }
}
else {
Channel.fromPath(params.fastqlist)
        .splitCsv()
        .map { reads ->
                [file(reads[0]), file(reads[1])]
        }.set { fastq }
}
process trimmomatic {
        publishDir "trimmed_fastq", pattern:'*.*', mode:'copy'
        input:
	set file(read1), file(read2) from fastq
        output:
        file "*_trim.fastq.gz" into trimmed_fastq
        script:
        if (params.reads == "SE") {
        """
	echo true
        """
	}
	else {
	"""
        trim1=\$(echo $read1 | sed 's/.fastq.gz/_trim.fastq.gz/g' | sed 's/.fq/_trim.fastq.gz/g')
        trim2=\$(echo $read2 | sed 's/.fastq.gz/_trim.fastq.gz/g' | sed 's/.fq/_trim.fastq.gz/g')
        sample=\$(echo $read2 | sed 's/_R1_001.fastq.gz//g' | sed 's/_R1.fastq.gz//g' | sed 's/_R1.fq//g')
        java -jar $path2trim PE -threads $task.cpus $read1 $read2 \$trim1 ${read1}_unpaired \$trim2 ${read2}_unpaired ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$readlength CROP:$readlength &> \${sample}_log.txt
        """
	}
}
process star {
        publishDir "star_out", pattern:'*.*' , mode: 'copy'
        input:
	set file(trim1), file(trim2) from trimmed_fastq
        output:
        file "*.bam" into star_out1, star_out2
        script:
        """
	module load STAR/2.5.3
        module load samtools/1.3.1
        sample=\$(echo $trim1 | sed 's/_R1_001_trim.fastq.gz//g' | sed 's/_R1_trim.fastq.gz//g' )
        overhang=\$(($readlength -1))
        STAR --genomeDir $genomeDir --readFilesIn $trim1 $trim2 --readFilesCommand zcat --readMatesLengthsIn NotEqual --outFileNamePrefix \$sample. --runThreadN $task.cpus -sjdbGTFfile $refgtf --sjdbOverhang \$overhang --alignSJoverhangMin 8 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax 2 --outFilterMultimapNmax 20  --outReadsUnmapped Fastx --alignMatesGapMax 1000000 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000  --outBAMsortingThreadN $task.cpus --outFilterType BySJout --twopassMode Basic --alignEndsType EndToEnd
        """
}
process bam2bw {
        publishDir "bw_files", pattern: '*fin.bw', mode: 'copy'
        input:
	file(sample) from star_out1
        file(chrLen) from genomeDirchr
        output:
        file "*fin.bw" into bw_files
        script:
        """
	module load compsci
        module load python/2.7.3
        module load samtools/1.3.1
        sample2=\$(basename $sample .Aligned.sortedByCoord.out.bam)
        samtools index ${sample}
        bamCoverage -b ${sample} -o \$sample2.bw
        /opt/compsci/kent/bin/bigWigToWig \$sample2.bw \$sample2.tmp1.wig
        grep -v bedGraph \$sample2.tmp1.wig > \$sample2.fin.wig
        /opt/compsci/kent/bin/wigToBigWig \$sample2.fin.wig ${chrLen} \$sample2.fin.bw
        """
}
process stringtie {
        publishDir "stringtie_out", pattern: '*.gtf', mode: 'copy'
        input:
	file(sample) from star_out2
        output:
        file('*.for.rMATS.gtf') into stringtie_out
        file('*.gtf')
        script:
        """
	module load compsci
        samplename=\$(basename $sample .Aligned.sortedByCoord.out.bam)
        /projects/anczukow-lab/tools/stringtie-1.3.6/stringtie $sample -G $refgtf -o \$samplename.for.DGE.gtf --rf -e -a 8 -p $task.cpus
        /projects/anczukow-lab/tools/stringtie-1.3.6/stringtie $sample -G $refgtf -o \$samplename.for.rMATS.gtf --rf -a 8 -p $task.cpus
        """
}
process stringtiemerge {
        publishDir 'merged_gtf', pattern:'*.gtf', mode:'copy'
        input:
	file refgtf
        file('*.gtf') from stringtie_out.collect()
        val outgtf
        output:
        file('*.gtf')
        """
	module load compsci
        module load it-modules
        ls -1  *.gtf >  list_of_gtf_files_to_merge.txt
        /projects/anczukow-lab/tools/stringtie-1.3.6/stringtie --merge -p 8 -o  $outgtf -G $refgtf list_of_gtf_files_to_merge.txt
        """
}