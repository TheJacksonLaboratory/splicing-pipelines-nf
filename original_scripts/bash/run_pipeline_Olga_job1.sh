#!/bin/sh
#PBS -q batch
#PBS -l nodes=1:ppn=14
#PBS -l walltime=72:00:00

# use script: qsub -d . -v reads="PE",stranded="yes",dataset="Dataset_3_SRSF1mutants_MCF10A",readlength="48",outdir="/projects/anczukow-lab/yuriem/test_pipeline",threads="14",stargenome="/projects/anczukow-lab/human_transcriptome_reference_files/star_overhang_47",gtf="/data/yuriem/bin/genomes/Annotation/Homo_sapiens.GRCh38.98.gtf",mismatch="2" run_pipeline_Olga_with_arguments.pbs

module load fastqc/0.11.5  
module load java
module load samtools/1.5
module load python/2.7.3
module load compsci


if [ "$reads" == "PE" ]; then
	read1=$(cat ${outdir}/list1.txt | head -n $PBS_ARRAYID | tail -n 1)

	read2=$(echo $read1 | sed 's/_R1.f/_R2.f/g' | sed 's/_R1_/_R2_/g')

	read1_new=$(echo $read1 | rev | cut -d / -f1 | rev )
	read2_new=$(echo $read2 | rev | cut -d / -f1 | rev )
	out_read1=$(echo $read1_new | sed 's/_R1_/_R1_unpaired_/g' )
	out_read2=$(echo $read2_new | sed 's/_R2_/_R2_unpaired_/g' )

	overhang="$(($readlength -1))"
	
	echo $read1 $read2
	echo $read1_new $read2_new
#FASTQC

	fastqc $read1  --casava -t $threads --outdir=${outdir}/QC_raw
	fastqc $read2 --casava -t $threads --outdir=${outdir}/QC_raw

#TRIMMOMATIC

	java -jar /projects/anczukow-lab/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $threads -phred33 ${read1} ${read2} ${outdir}/trimmed/${read1_new} ${outdir}/trimmed/${out_read1} ${outdir}/trimmed/${read2_new} ${outdir}/trimmed/${out_read2} ILLUMINACLIP:/projects/anczukow-lab/tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${readlength} CROP:${readlength}

#STAR

	sample=$(echo $read1_new | sed 's/_R1.fastq.gz//g' | sed 's/_R1_001.fastq.gz//g')

	if [ $stranded == "yes" ]; then
		/projects/anczukow-lab/tools/STAR-2.7.3a/bin/Linux_x86_64_static/STAR --genomeDir $stargenome --readFilesIn ${outdir}/trimmed/${read1_new} ${outdir}/trimmed/${read2_new} --readMatesLengthsIn NotEqual --outFileNamePrefix ${outdir}/star_mapped/${sample}  --runThreadN $threads --readFilesCommand zcat --sjdbGTFfile ${gtf} --sjdbOverhang ${overhang} --alignSJoverhangMin 8 --outFilterMismatchNmax $mismatch --outFilterMultimapNmax 20 --alignMatesGapMax 1000000 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000  --outBAMsortingThreadN $threads --outFilterType BySJout --twopassMode Basic --alignEndsType EndToEnd --outWigType wiggle

# make bigWig files

#		wigToBigWig ${outdir}/star_mapped/${sample}*.wig $chromfile ${outdir}/star_mapped/${sample}.bw

		samtools index ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam

		bamCoverage -b ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -o ${outdir}/star_mapped/${sample}_old.bw  


#Stringtie

		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}.gtf --rf -a 8 -p $threads

#Stringtie DEG
		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}_for_DGE.gtf --rf -a 8 -e -p $threads

#Not stranded
	else
		/projects/anczukow-lab/tools/STAR-2.7.3a/bin/Linux_x86_64_static/STAR --genomeDir $stargenome --readFilesIn ${outdir}/trimmed/${read1_new} ${outdir}/trimmed/${read2_new} --readMatesLengthsIn NotEqual --outFileNamePrefix ${outdir}/star_mapped/${sample}  --runThreadN $threads --readFilesCommand zcat --sjdbGTFfile ${gtf} --sjdbOverhang ${overhang} --alignSJoverhangMin 8 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax $mismatch --outFilterMultimapNmax 20 --alignMatesGapMax 1000000 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000  --outBAMsortingThreadN $threads --outFilterType BySJout --twopassMode Basic --alignEndsType EndToEnd

#Stringtie

		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}.gtf -a 8 -p $threads

#Stringtie DEG
		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}_for_DGE.gtf -a 8 -e -p $threads

	fi

#Single end

else

	read1=$(cat ${outdir}/list1.txt | head -n $PBS_ARRAYID | tail -n 1)

	read1_new=$(echo $read1 | rev | cut -d / -f1 | rev)

	overhang="$(($readlength -1))"

#FASTQC

	fastqc $read1  --casava -t $threads --outdir=${outdir}/QC_raw

#TRIMMOMATIC
	
	java -jar /projects/anczukow-lab/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $threads -phred33 ${read1} ${outdir}/trimmed/${read1_new} ILLUMINACLIP:/projects/anczukow-lab/tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${readlength} CROP:${readlength}

#STAR

	sample=$(echo $read1_new | sed 's/_R1_001.fastq.gz//g')

#Stranded
	if [ $stranded == "yes" ]; then
#STAR
		/projects/anczukow-lab/tools/STAR-2.7.3a/bin/Linux_x86_64_static/STAR --genomeDir $stargenome --readFilesIn ${outdir}/trimmed/${read1_new} --readMatesLengthsIn NotEqual --outFileNamePrefix ${outdir}/star_mapped/${sample}  --runThreadN $threads --readFilesCommand zcat  --sjdbGTFfile ${gtf} --sjdbOverhang ${overhang} --alignSJoverhangMin 8 --outFilterMismatchNmax $mismatch --outFilterMultimapNmax 20 --alignMatesGapMax 1000000 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000  --outBAMsortingThreadN $threads --outFilterType BySJout --twopassMode Basic --alignEndsType EndToEnd


#Stringtie

		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}.gtf --rf -a 8 -p $threads

#Stringtie for DEG
		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}_for_DGE.gtf --rf -a 8 -e -p $threads

#Not stranded
	else
#STAR

		/projects/anczukow-lab/tools/STAR-2.7.3a/bin/Linux_x86_64_static/STAR --genomeDir $stargenome --readFilesIn ${outdir}/trimmed/${read1_new} --readMatesLengthsIn NotEqual --outFileNamePrefix ${outdir}/star_mapped/${sample}  --runThreadN $threads --readFilesCommand zcat --sjdbGTFfile ${gtf} --sjdbOverhang ${overhang} --alignSJoverhangMin 8 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax $mismatch --outFilterMultimapNmax 20 --alignMatesGapMax 1000000 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000  --outBAMsortingThreadN $threads --outFilterType BySJout --twopassMode Basic --alignEndsType EndToEnd

#Stringtie
		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}.gtf -a 8 -p $threads

#Stringtie for DEG
		/projects/anczukow-lab/tools/stringtie-2.0.6/stringtie ${outdir}/star_mapped/${sample}Aligned.sortedByCoord.out.bam -G $gtf -o ${outdir}/star_mapped/${sample}_for_DGE.gtf -a 8 -e -p $threads

	fi	
fi

echo "Part1 is complete!"
