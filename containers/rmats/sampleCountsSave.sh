#!/bin/bash
#
# 2017 October 16
# 2019 July 5 - ADM - updated to work with rMATS 4.0.2 
#
# Purpose:  To save the temporary files that are created running RNASeq-MATS.py
#           Saving these 10 files will enable the creation of the final matrix
#           so that subsequent analysis can be done.
#
#           For each GTF, TCGA nodup.bam files that had been created
#           according to the following method.  Each RNASeq fastq file is
#           1.  Trimmed to the lowest common denominator length (48 BP) and
#           2.  Aligned using hisat2 as of this date to GRCh38 genome assembly,
#           3.  duplicates removed
#           4.  and coordinate sorted.
#           These files are heretoafter refered to as our nodup.bam files.
#
#           Using this resource -- which is constructed in a previous step, these
#           files are compared to the provided GTF.
#
#           rMATS creates 5 files that we use which are the junctions that define
#           5 specific kinds of alternative splicing events.
#
#           1. Alternative 3' Splice Site
#              -- the junction IDs, coordinates and involved exons provided in fromGTF.A3SS.txt
#           2. Alternative 5' Splice Site
#              -- the junction IDs, coordinates and involved exons provided in fromGTF.A5SS.txt
#           3. Mutually exclusive exons
#              -- the junction IDS, coordinates and involved exons provided in fromGTF.MXE.txt
#           4. Retention Introns
#              -- the junction IDs, coordinates and involved exons provided in fromGTF.RI.txt
#           5. Skipped Exons
#              -- the junction IDs, coordiantes and involved exons provided in fromGTF.SE.txt
#
#           Note that the columns 6 and 7 are the normalized length for inclusion and skipped exons.
#           The field 7 will always be trim length -1 == in our case the least common denominator is 48
#           So the trim length is 48 bp -- so field 7 is 47 for exon skipping events -- more variable for
#           mutually exclusive exons as there will always be a variable length for the included exon.
#           In any case they are never zero.
#
# Algorithm Design
#
# For each rMATS comparison -- which was made at the single comparison level, there is a repeated directory
# ASEvents -- in this directory are the following 5 tab delimited filenames together with their header indicating content
#
# 2019 July 5 - ADM - update.  There is no longer an ASEvents directory or a temp file directory
#             The output of the ASEvents is kept in the same directory as the results
#
# One copy of this will be preserved
#
#   fromGTF.A3SS.txt (tabs replaced by space)
#   ID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
#
#   fromGTF.A5SS.txt
#   ID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
#
#   fromGTF.MXE.txt
#   ID GeneID geneSymbol chr strand 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
#
#   fromGTF.RI.txt
#   ID GeneID geneSymbol chr strand riExonStart_0base riExonEnd upstreamES upstreamEE downstreamES downstreamEE
#
#   fromGTF.SE.txt
#   ID  GeneID geneSymbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
#
# They contain the junction id's for each of these events -- these are the same for each GTF.
#
#
# Anne Deslattes Mays, PhD
#
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# 2017 October 18 -- this worked in our test but the summary.txt file is not being found -- here is hoping this works
# the summary.txt file is at the highest level directory
#   -- we can extract from this file the name of the two samples
# sample1=$(sed '3q;d' $dir"/summary.txt" | awk -F "/" '{print $5}' | awk -F "." '{print $1}')
# sample2=$(sed '4q;d' $dir"/summary.txt" | awk -F "/" '{print $5}' | awk -F "." '{print $1}')
#

# 2019 July 5 -- there is no longer a summary.txt file or a temp directory
#     Sample names can be read from the command line -- fed from the rMATS 4.0.2 script rnaseq_mats_sample1_sample2_tar_4.0.2.py
#
#
#------------------------------------------------------------------------------------------------------
dirraw=$1
dir=$1"/"
sample1=$2
sample2=$3
echo "output directory = "
echo $dir

echo "sample1 ="
echo $sample1

echo "sample2 ="
echo $sample2
#------------------------------------------------------------------------------------------------------
#
# output files names (endings)
#
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#
# A3SS
#
#------------------------------------------------------------------------------------------------------
a3ss_jc_idend=".a3ss.jc.id.txt"
a3ss_jc_incend=".a3ss.jc.inc.txt"
a3ss_jc_ijcend=".a3ss.jc.ijc.txt"
a3ss_jc_sjcend=".a3ss.jc.sjc.txt"
a3ss_jc_inclenend=".a3ss.jc.inclen.txt"
a3ss_jc_skiplenend=".a3ss.jc.skiplen.txt"

a3ss_jcec_idend=".a3ss.jcec.id.txt"
a3ss_jcec_incend=".a3ss.jcec.inc.txt"
a3ss_jcec_ijcend=".a3ss.jcec.ijc.txt"
a3ss_jcec_sjcend=".a3ss.jcec.sjc.txt"
a3ss_jcec_inclenend=".a3ss.jcec.inclen.txt"
a3ss_jcec_skiplenend=".a3ss.jcec.skiplen.txt"
#------------------------------------------------------------------------------------------------------
#
# A5SS
#
#------------------------------------------------------------------------------------------------------
a5ss_jc_idend=".a5ss.jc.id.txt"
a5ss_jc_incend=".a5ss.jc.inc.txt"
a5ss_jc_ijcend=".a5ss.jc.ijc.txt"
a5ss_jc_sjcend=".a5ss.jc.sjc.txt"
a5ss_jc_inclenend=".a5ss.jc.inclen.txt"
a5ss_jc_skiplenend=".a5ss.jc.skiplen.txt"

a5ss_jcec_idend=".a5ss.jcec.id.txt"
a5ss_jcec_incend=".a5ss.jcec.inc.txt"
a5ss_jcec_ijcend=".a5ss.jcec.ijc.txt"
a5ss_jcec_sjcend=".a5ss.jcec.sjc.txt"
a5ss_jcec_inclenend=".a5ss.jcec.inclen.txt"
a5ss_jcec_skiplenend=".a5ss.jcec.skiplen.txt"
#------------------------------------------------------------------------------------------------------
#
# MXE
#
#------------------------------------------------------------------------------------------------------
mxe_jc_idend=".mxe.jc.id.txt"
mxe_jc_incend=".mxe.jc.inc.txt"
mxe_jc_ijcend=".mxe.jc.ijc.txt"
mxe_jc_sjcend=".mxe.jc.sjc.txt"
mxe_jc_inclenend=".mxe.jc.inclen.txt"
mxe_jc_skiplenend=".mxe.jc.skiplen.txt"

mxe_jcec_idend=".mxe.jcec.id.txt"
mxe_jcec_incend=".mxe.jcec.inc.txt"
mxe_jcec_ijcend=".mxe.jcec.ijc.txt"
mxe_jcec_sjcend=".mxe.jcec.sjc.txt"
mxe_jcec_inclenend=".mxe.jcec.inclen.txt"
mxe_jcec_skiplenend=".mxe.jcec.skiplen.txt"
#------------------------------------------------------------------------------------------------------
#
# RI
#
#------------------------------------------------------------------------------------------------------
ri_jc_idend=".ri.jc.id.txt"
ri_jc_incend=".ri.jc.inc.txt"
ri_jc_ijcend=".ri.jc.ijc.txt"
ri_jc_sjcend=".ri.jc.sjc.txt"
ri_jc_inclenend=".ri.jc.inclen.txt"
ri_jc_skiplenend=".ri.jc.skiplen.txt"

ri_jcec_idend=".ri.jcec.id.txt"
ri_jcec_incend=".ri.jcec.inc.txt"
ri_jcec_ijcend=".ri.jcec.ijc.txt"
ri_jcec_sjcend=".ri.jcec.sjc.txt"
ri_jcec_inclenend=".ri.jcec.inclen.txt"
ri_jcec_skiplenend=".ri.jcec.skiplen.txt"
#------------------------------------------------------------------------------------------------------
#
# SE
#
#------------------------------------------------------------------------------------------------------
se_jc_idend=".se.jc.id.txt"
se_jc_incend=".se.jc.inc.txt"
se_jc_ijcend=".se.jc.ijc.txt"
se_jc_sjcend=".se.jc.sjc.txt"
se_jc_inclenend=".se.jc.inclen.txt"
se_jc_skiplenend=".se.jc.skiplen.txt"

se_jcec_idend=".se.jcec.id.txt"
se_jcec_incend=".se.jcec.inc.txt"
se_jcec_ijcend=".se.jcec.ijc.txt"
se_jcec_sjcend=".se.jcec.sjc.txt"
se_jcec_inclenend=".se.jcec.inclen.txt"
se_jcec_skiplenend=".se.jcec.skiplen.txt"
#------------------------------------------------------------------------------------------------------
#
#  sample1 - file names
#
#------------------------------------------------------------------------------------------------------
#
# A3SS
#
#------------------------------------------------------------------------------------------------------
a3ss_jc_idend_sample1=$dir$sample1$a3ss_jc_idend
echo "a3ss_jc_idend_sample1="
echo $a3ss_jc_idend_sample1
a3ss_jc_incend_sample1=$dir$sample1$a3ss_jc_incend
echo "a3ss_jc_incend_sample1="
echo $a3ss_jc_incend_sample1
a3ss_jc_ijcend_sample1=$dir$sample1$a3ss_jc_ijcend
echo "a3ss_jc_ijcend_sample1="
echo $a3ss_jc_ijcend_sample1
a3ss_jc_sjcend_sample1=$dir$sample1$a3ss_jc_sjcend
echo "a3ss_jc_sjcend_sample1="
echo $a3ss_jc_sjcend_sample1
a3ss_jc_inclenend_sample1=$dir$sample1$a3ss_jc_inclenend
echo "a3ss_jc_inclenend_sample1="
echo $a3ss_jc_inclenend_sample1
a3ss_jc_skiplenend_sample1=$dir$sample1$a3ss_jc_skiplenend
echo "a3ss_jc_skiplenend_sample1="
echo $a3ss_jc_skiplenend_sample1

a3ss_jcec_idend_sample1=$dir$sample1$a3ss_jcec_idend
echo "a3ss_jcec_idend_sample1="
echo $a3ss_jcec_idend_sample1
a3ss_jcec_incend_sample1=$dir$sample1$a3ss_jcec_incend
echo "a3ss_jcec_incend_sample1="
echo $a3ss_jcec_incend_sample1
a3ss_jcec_ijcend_sample1=$dir$sample1$a3ss_jcec_ijcend
echo "a3ss_jcec_ijcend_sample1="
echo $a3ss_jcec_ijcend_sample1
a3ss_jcec_sjcend_sample1=$dir$sample1$a3ss_jcec_sjcend
echo "a3ss_jcec_sjcend_sample1="
echo $a3ss_jcec_sjcend_sample1
a3ss_jcec_inclenend_sample1=$dir$sample1$a3ss_jcec_inclenend
echo "a3ss_jcec_inclenend_sample1="
echo $a3ss_jcec_inclenend_sample1
a3ss_jcec_skiplenend_sample1=$dir$sample1$a3ss_jcec_skiplenend
echo "a3ss_jcec_skiplenend_sample1="
echo $a3ss_jcec_skiplenend_sample1
#------------------------------------------------------------------------------------------------------
#
# A5SS
#
#------------------------------------------------------------------------------------------------------
a5ss_jc_idend_sample1=$dir$sample1$a5ss_jc_idend
echo "a5ss_jc_idend_sample1="
echo $a5ss_jc_idend_sample1
a5ss_jc_incend_sample1=$dir$sample1$a5ss_jc_incend
echo "a5ss_jc_incend_sample1="
echo $a5ss_jc_incend_sample1
a5ss_jc_ijcend_sample1=$dir$sample1$a5ss_jc_ijcend
echo "a5ss_jc_ijcend_sample1="
echo $a5ss_jc_ijcend_sample1
a5ss_jc_sjcend_sample1=$dir$sample1$a5ss_jc_sjcend
echo "a5ss_jc_sjcend_sample1="
echo $a5ss_jc_sjcend_sample1
a5ss_jc_inclenend_sample1=$dir$sample1$a5ss_jc_inclenend
echo "a5ss_jc_inclenend_sample1="
echo $a5ss_jc_inclenend_sample1
a5ss_jc_skiplenend_sample1=$dir$sample1$a5ss_jc_skiplenend
echo "a5ss_jc_skiplenend_sample1="
echo $a5ss_jc_skiplenend_sample1

a5ss_jcec_idend_sample1=$dir$sample1$a5ss_jcec_idend
echo "a5ss_jcec_idend_sample1="
echo $a5ss_jcec_idend_sample1
a5ss_jcec_incend_sample1=$dir$sample1$a5ss_jcec_incend
echo "a5ss_jcec_incend_sample1="
echo $a5ss_jcec_incend_sample1
a5ss_jcec_ijcend_sample1=$dir$sample1$a5ss_jcec_ijcend
echo "a5ss_jcec_ijcend_sample1="
echo $a5ss_jcec_ijcend_sample1
a5ss_jcec_sjcend_sample1=$dir$sample1$a5ss_jcec_sjcend
echo "a5ss_jcec_sjcend_sample1="
echo $a5ss_jcec_sjcend_sample1
a5ss_jcec_inclenend_sample1=$dir$sample1$a5ss_jcec_inclenend
echo "a5ss_jcec_inclenend_sample1="
echo $a5ss_jcec_inclenend_sample1
a5ss_jcec_skiplenend_sample1=$dir$sample1$a5ss_jcec_skiplenend
echo "a5ss_jcec_skiplenend_sample1="
echo $a5ss_jcec_skiplenend_sample1
#------------------------------------------------------------------------------------------------------
#
# MXE
#
#------------------------------------------------------------------------------------------------------
mxe_jc_idend_sample1=$dir$sample1$mxe_jc_idend
echo "mxe_jc_idend_sample1="
echo $mxe_jc_idend_sample1
mxe_jc_incend_sample1=$dir$sample1$mxe_jc_incend
echo "mxe_jc_incend_sample1="
echo $mxe_jc_incend_sample1
mxe_jc_ijcend_sample1=$dir$sample1$mxe_jc_ijcend
echo "mxe_jc_ijcend_sample1="
echo $mxe_jc_ijcend_sample1
mxe_jc_sjcend_sample1=$dir$sample1$mxe_jc_sjcend
echo "mxe_jc_sjcend_sample1="
echo $mxe_jc_sjcend_sample1
mxe_jc_inclenend_sample1=$dir$sample1$mxe_jc_inclenend
echo "mxe_jc_inclenend_sample1="
echo $mxe_jc_inclenend_sample1
mxe_jc_skiplenend_sample1=$dir$sample1$mxe_jc_skiplenend
echo "mxe_jc_skiplenend_sample1="
echo $mxe_jc_skiplenend_sample1

mxe_jcec_idend_sample1=$dir$sample1$mxe_jcec_idend
echo "mxe_jcec_idend_sample1="
echo $mxe_jcec_idend_sample1
mxe_jcec_incend_sample1=$dir$sample1$mxe_jcec_incend
echo "mxe_jcec_incend_sample1="
echo $mxe_jcec_incend_sample1
mxe_jcec_ijcend_sample1=$dir$sample1$mxe_jcec_ijcend
echo "mxe_jcec_ijcend_sample1="
echo $mxe_jcec_ijcend_sample1
mxe_jcec_sjcend_sample1=$dir$sample1$mxe_jcec_sjcend
echo "mxe_jcec_sjcend_sample1="
echo $mxe_jcec_sjcend_sample1
mxe_jcec_inclenend_sample1=$dir$sample1$mxe_jcec_inclenend
echo "mxe_jcec_inclenend_sample1="
echo $mxe_jcec_inclenend_sample1
mxe_jcec_skiplenend_sample1=$dir$sample1$mxe_jcec_skiplenend
echo "mxe_jcec_skiplenend_sample1="
echo $mxe_jcec_skiplenend_sample1
#------------------------------------------------------------------------------------------------------
#
# RI
#
#------------------------------------------------------------------------------------------------------
ri_jc_idend_sample1=$dir$sample1$ri_jc_idend
echo "ri_jc_idend_sample1="
echo $ri_jc_idend_sample1
ri_jc_incend_sample1=$dir$sample1$ri_jc_incend
echo "ri_jc_incend_sample1="
echo $ri_jc_incend_sample1
ri_jc_ijcend_sample1=$dir$sample1$ri_jc_ijcend
echo "ri_jc_ijcend_sample1="
echo $ri_jc_ijcend_sample1
ri_jc_sjcend_sample1=$dir$sample1$ri_jc_sjcend
echo "ri_jc_sjcend_sample1="
echo $ri_jc_sjcend_sample1
ri_jc_inclenend_sample1=$dir$sample1$ri_jc_inclenend
echo "ri_jc_inclenend_sample1="
echo $ri_jc_inclenend_sample1
ri_jc_skiplenend_sample1=$dir$sample1$ri_jc_skiplenend
echo "ri_jc_skiplenend_sample1="
echo $ri_jc_skiplenend_sample1

ri_jcec_idend_sample1=$dir$sample1$ri_jcec_idend
echo "ri_jcec_idend_sample1="
echo $ri_jcec_idend_sample1
ri_jcec_incend_sample1=$dir$sample1$ri_jcec_incend
echo "ri_jcec_incend_sample1="
echo $ri_jcec_incend_sample1
ri_jcec_ijcend_sample1=$dir$sample1$ri_jcec_ijcend
echo "ri_jcec_ijcend_sample1="
echo $ri_jcec_ijcend_sample1
ri_jcec_sjcend_sample1=$dir$sample1$ri_jcec_sjcend
echo "ri_jcec_sjcend_sample1="
echo $ri_jcec_sjcend_sample1
ri_jcec_inclenend_sample1=$dir$sample1$ri_jcec_inclenend
echo "ri_jcec_inclenend_sample1="
echo $ri_jcec_inclenend_sample1
ri_jcec_skiplenend_sample1=$dir$sample1$ri_jcec_skiplenend
echo "ri_jcec_skiplenend_sample1="
echo $ri_jcec_skiplenend_sample1
#------------------------------------------------------------------------------------------------------
#
# SE
#
#------------------------------------------------------------------------------------------------------
se_jc_idend_sample1=$dir$sample1$se_jc_idend
echo "se_jc_idend_sample1="
echo $se_jc_idend_sample1
se_jc_incend_sample1=$dir$sample1$se_jc_incend
echo "se_jc_incend_sample1="
echo $se_jc_incend_sample1
se_jc_ijcend_sample1=$dir$sample1$se_jc_ijcend
echo "se_jc_ijcend_sample1="
echo $se_jc_ijcend_sample1
se_jc_sjcend_sample1=$dir$sample1$se_jc_sjcend
echo "se_jc_sjcend_sample1="
echo $se_jc_sjcend_sample1
se_jc_inclenend_sample1=$dir$sample1$se_jc_inclenend
echo "se_jc_inclenend_sample1="
echo $se_jc_inclenend_sample1
se_jc_skiplenend_sample1=$dir$sample1$se_jc_skiplenend
echo "se_jc_skiplenend_sample1="
echo $se_jc_skiplenend_sample1

se_jcec_idend_sample1=$dir$sample1$se_jcec_idend
echo "se_jcec_idend_sample1="
echo $se_jcec_idend_sample1
se_jcec_incend_sample1=$dir$sample1$se_jcec_incend
echo "se_jcec_incend_sample1="
echo $se_jcec_incend_sample1
se_jcec_ijcend_sample1=$dir$sample1$se_jcec_ijcend
echo "se_jcec_ijcend_sample1="
echo $se_jcec_ijcend_sample1
se_jcec_sjcend_sample1=$dir$sample1$se_jcec_sjcend
echo "se_jcec_sjcend_sample1="
echo $se_jcec_sjcend_sample1
se_jcec_inclenend_sample1=$dir$sample1$se_jcec_inclenend
echo "se_jcec_inclenend_sample1="
echo $se_jcec_inclenend_sample1
se_jcec_skiplenend_sample1=$dir$sample1$se_jcec_skiplenend
echo "se_jcec_skiplenend_sample1="
echo $se_jcec_skiplenend_sample1
#------------------------------------------------------------------------------------------------------
#
#  sample2 - file names
#
#------------------------------------------------------------------------------------------------------
#
# A3SS
#
#------------------------------------------------------------------------------------------------------
a3ss_jc_idend_sample2=$dir$sample2$a3ss_jc_idend
echo "a3ss_jc_idend_sample2="
echo $a3ss_jc_idend_sample2
a3ss_jc_incend_sample2=$dir$sample2$a3ss_jc_incend
echo "a3ss_jc_incend_sample2="
echo $a3ss_jc_incend_sample2
a3ss_jc_ijcend_sample2=$dir$sample2$a3ss_jc_ijcend
echo "a3ss_jc_ijcend_sample2="
echo $a3ss_jc_ijcend_sample2
a3ss_jc_sjcend_sample2=$dir$sample2$a3ss_jc_sjcend
echo "a3ss_jc_sjcend_sample2="
echo $a3ss_jc_sjcend_sample2
a3ss_jc_inclenend_sample2=$dir$sample2$a3ss_jc_inclenend
echo "a3ss_jc_inclenend_sample2="
echo $a3ss_jc_inclenend_sample2
a3ss_jc_skiplenend_sample2=$dir$sample2$a3ss_jc_skiplenend
echo "a3ss_jc_skiplenend_sample2="
echo $a3ss_jc_skiplenend_sample2

a3ss_jcec_idend_sample2=$dir$sample2$a3ss_jcec_idend
echo "a3ss_jcec_idend_sample2="
echo $a3ss_jcec_idend_sample2
a3ss_jcec_incend_sample2=$dir$sample2$a3ss_jcec_incend
echo "a3ss_jcec_incend_sample2="
echo $a3ss_jcec_incend_sample2
a3ss_jcec_ijcend_sample2=$dir$sample2$a3ss_jcec_ijcend
echo "a3ss_jcec_ijcend_sample2="
echo $a3ss_jcec_ijcend_sample2
a3ss_jcec_sjcend_sample2=$dir$sample2$a3ss_jcec_sjcend
echo "a3ss_jcec_sjcend_sample2="
echo $a3ss_jcec_sjcend_sample2
a3ss_jcec_inclenend_sample2=$dir$sample2$a3ss_jcec_inclenend
echo "a3ss_jcec_inclenend_sample2="
echo $a3ss_jcec_inclenend_sample2
a3ss_jcec_skiplenend_sample2=$dir$sample2$a3ss_jcec_skiplenend
echo "a3ss_jcec_skiplenend_sample2="
echo $a3ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# A5SS
#
#------------------------------------------------------------------------------------------------------
a5ss_jc_idend_sample2=$dir$sample2$a5ss_jc_idend
echo "a5ss_jc_idend_sample2="
echo $a5ss_jc_idend_sample2
a5ss_jc_incend_sample2=$dir$sample2$a5ss_jc_incend
echo "a5ss_jc_incend_sample2="
echo $a5ss_jc_incend_sample2
a5ss_jc_ijcend_sample2=$dir$sample2$a5ss_jc_ijcend
echo "a5ss_jc_ijcend_sample2="
echo $a5ss_jc_ijcend_sample2
a5ss_jc_sjcend_sample2=$dir$sample2$a5ss_jc_sjcend
echo "a5ss_jc_sjcend_sample2="
echo $a5ss_jc_sjcend_sample2
a5ss_jc_inclenend_sample2=$dir$sample2$a5ss_jc_inclenend
echo "a5ss_jc_inclenend_sample2="
echo $a5ss_jc_inclenend_sample2
a5ss_jc_skiplenend_sample2=$dir$sample2$a5ss_jc_skiplenend
echo "a5ss_jc_skiplenend_sample2="
echo $a5ss_jc_skiplenend_sample2

a5ss_jcec_idend_sample2=$dir$sample2$a5ss_jcec_idend
echo "a5ss_jcec_idend_sample2="
echo $a5ss_jcec_idend_sample2
a5ss_jcec_incend_sample2=$dir$sample2$a5ss_jcec_incend
echo "a5ss_jcec_incend_sample2="
echo $a5ss_jcec_incend_sample2
a5ss_jcec_ijcend_sample2=$dir$sample2$a5ss_jcec_ijcend
echo "a5ss_jcec_ijcend_sample2="
echo $a5ss_jcec_ijcend_sample2
a5ss_jcec_sjcend_sample2=$dir$sample2$a5ss_jcec_sjcend
echo "a5ss_jcec_sjcend_sample2="
echo $a5ss_jcec_sjcend_sample2
a5ss_jcec_inclenend_sample2=$dir$sample2$a5ss_jcec_inclenend
echo "a5ss_jcec_inclenend_sample2="
echo $a5ss_jcec_inclenend_sample2
a5ss_jcec_skiplenend_sample2=$dir$sample2$a5ss_jcec_skiplenend
echo "a5ss_jcec_skiplenend_sample2="
echo $a5ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# MXE
#
#------------------------------------------------------------------------------------------------------
mxe_jc_idend_sample2=$dir$sample2$mxe_jc_idend
echo "mxe_jc_idend_sample2="
echo $mxe_jc_idend_sample2
mxe_jc_incend_sample2=$dir$sample2$mxe_jc_incend
echo "mxe_jc_incend_sample2="
echo $mxe_jc_incend_sample2
mxe_jc_ijcend_sample2=$dir$sample2$mxe_jc_ijcend
echo "mxe_jc_ijcend_sample2="
echo $mxe_jc_ijcend_sample2
mxe_jc_sjcend_sample2=$dir$sample2$mxe_jc_sjcend
echo "mxe_jc_sjcend_sample2="
echo $mxe_jc_sjcend_sample2
mxe_jc_inclenend_sample2=$dir$sample2$mxe_jc_inclenend
echo "mxe_jc_inclenend_sample2="
echo $mxe_jc_inclenend_sample2
mxe_jc_skiplenend_sample2=$dir$sample2$mxe_jc_skiplenend
echo "mxe_jc_skiplenend_sample2="
echo $mxe_jc_skiplenend_sample2

mxe_jcec_idend_sample2=$dir$sample2$mxe_jcec_idend
echo "mxe_jcec_idend_sample2="
echo $mxe_jcec_idend_sample2
mxe_jcec_incend_sample2=$dir$sample2$mxe_jcec_incend
echo "mxe_jcec_incend_sample2="
echo $mxe_jcec_incend_sample2
mxe_jcec_ijcend_sample2=$dir$sample2$mxe_jcec_ijcend
echo "mxe_jcec_ijcend_sample2="
echo $mxe_jcec_ijcend_sample2
mxe_jcec_sjcend_sample2=$dir$sample2$mxe_jcec_sjcend
echo "mxe_jcec_sjcend_sample2="
echo $mxe_jcec_sjcend_sample2
mxe_jcec_inclenend_sample2=$dir$sample2$mxe_jcec_inclenend
echo "mxe_jcec_inclenend_sample2="
echo $mxe_jcec_inclenend_sample2
mxe_jcec_skiplenend_sample2=$dir$sample2$mxe_jcec_skiplenend
echo "mxe_jcec_skiplenend_sample2="
echo $mxe_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# RI
#
#------------------------------------------------------------------------------------------------------
ri_jc_idend_sample2=$dir$sample2$ri_jc_idend
echo "ri_jc_idend_sample2="
echo $ri_jc_idend_sample2
ri_jc_incend_sample2=$dir$sample2$ri_jc_incend
echo "ri_jc_incend_sample2="
echo $ri_jc_incend_sample2
ri_jc_ijcend_sample2=$dir$sample2$ri_jc_ijcend
echo "ri_jc_ijcend_sample2="
echo $ri_jc_ijcend_sample2
ri_jc_sjcend_sample2=$dir$sample2$ri_jc_sjcend
echo "ri_jc_sjcend_sample2="
echo $ri_jc_sjcend_sample2
ri_jc_inclenend_sample2=$dir$sample2$ri_jc_inclenend
echo "ri_jc_inclenend_sample2="
echo $ri_jc_inclenend_sample2
ri_jc_skiplenend_sample2=$dir$sample2$ri_jc_skiplenend
echo "ri_jc_skiplenend_sample2="
echo $ri_jc_skiplenend_sample2

ri_jcec_idend_sample2=$dir$sample2$ri_jcec_idend
echo "ri_jcec_idend_sample2="
echo $ri_jcec_idend_sample2
ri_jcec_incend_sample2=$dir$sample2$ri_jcec_incend
echo "ri_jcec_incend_sample2="
echo $ri_jcec_incend_sample2
ri_jcec_ijcend_sample2=$dir$sample2$ri_jcec_ijcend
echo "ri_jcec_ijcend_sample2="
echo $ri_jcec_ijcend_sample2
ri_jcec_sjcend_sample2=$dir$sample2$ri_jcec_sjcend
echo "ri_jcec_sjcend_sample2="
echo $ri_jcec_sjcend_sample2
ri_jcec_inclenend_sample2=$dir$sample2$ri_jcec_inclenend
echo "ri_jcec_inclenend_sample2="
echo $ri_jcec_inclenend_sample2
ri_jcec_skiplenend_sample2=$dir$sample2$ri_jcec_skiplenend
echo "ri_jcec_skiplenend_sample2="
echo $ri_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# SE
#
#------------------------------------------------------------------------------------------------------
se_jc_idend_sample2=$dir$sample2$se_jc_idend
echo "se_jc_idend_sample2="
echo $se_jc_idend_sample2
se_jc_incend_sample2=$dir$sample2$se_jc_incend
echo "se_jc_incend_sample2="
echo $se_jc_incend_sample2
se_jc_ijcend_sample2=$dir$sample2$se_jc_ijcend
echo "se_jc_ijcend_sample2="
echo $se_jc_ijcend_sample2
se_jc_sjcend_sample2=$dir$sample2$se_jc_sjcend
echo "se_jc_sjcend_sample2="
echo $se_jc_sjcend_sample2
se_jc_inclenend_sample2=$dir$sample2$se_jc_inclenend
echo "se_jc_inclenend_sample2="
echo $se_jc_inclenend_sample2
se_jc_skiplenend_sample2=$dir$sample2$se_jc_skiplenend
echo "se_jc_skiplenend_sample2="
echo $se_jc_skiplenend_sample2

se_jcec_idend_sample2=$dir$sample2$se_jcec_idend
echo "se_jcec_idend_sample2="
echo $se_jcec_idend_sample2
se_jcec_incend_sample2=$dir$sample2$se_jcec_incend
echo "se_jcec_incend_sample2="
echo $se_jcec_incend_sample2
se_jcec_ijcend_sample2=$dir$sample2$se_jcec_ijcend
echo "se_jcec_ijcend_sample2="
echo $se_jcec_ijcend_sample2
se_jcec_sjcend_sample2=$dir$sample2$se_jcec_sjcend
echo "se_jcec_sjcend_sample2="
echo $se_jcec_sjcend_sample2
se_jcec_inclenend_sample2=$dir$sample2$se_jcec_inclenend
echo "se_jcec_inclenend_sample2="
echo $se_jcec_inclenend_sample2
se_jcec_skiplenend_sample2=$dir$sample2$se_jcec_skiplenend
echo "se_jcec_skiplenend_sample2="
echo $se_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# create sample specific matrices
#
#------------------------------------------------------------------------------------------------------
JC_RNASeq_A3SS_MATS=$dir"JC.raw.input.A3SS.txt"
echo "JC_RNASeq_A3SS_MATS="
echo $JC_RNASeq_A3SS_MATS
JC_RNASeq_A5SS_MATS=$dir"JC.raw.input.A5SS.txt"
echo "JC_RNASeq_A5SS_MATS="
echo $JC_RNASeq_A5SS_MATS
JC_RNASeq_MXE_MATS=$dir"JC.raw.input.MXE.txt"
echo "JC_RNASeq_MXE_MATS="
echo $JC_RNASeq_MXE_MATS
JC_RNASeq_RI_MATS=$dir"JC.raw.input.RI.txt"
echo "JC_RNASeq_RI_MATS="
echo $JC_RNASeq_RI_MATS
JC_RNASeq_SE_MATS=$dir"JC.raw.input.SE.txt"
echo "JC_RNASeq_SE_MATS="
echo $JC_RNASeq_SE_MATS

JCEC_RNASeq_A3SS_MATS=$dir"JCEC.raw.input.A3SS.txt"
echo "JCEC_RNASeq_A3SS_MATS="
echo $JCEC_RNASeq_A3SS_MATS
JCEC_RNASeq_A5SS_MATS=$dir"JCEC.raw.input.A5SS.txt"
echo "JCEC_RNASeq_A5SS_MATS="
echo $JCEC_RNASeq_A5SS_MATS
JCEC_RNASeq_MXE_MATS=$dir"JCEC.raw.input.MXE.txt"
echo "JCEC_RNASeq_MXE_MATS="
echo $JCEC_RNASeq_MXE_MATS
JCEC_RNASeq_RI_MATS=$dir"JCEC.raw.input.RI.txt"
echo "JCEC_RNASeq_RI_MATS="
echo $JCEC_RNASeq_RI_MATS
JCEC_RNASeq_SE_MATS=$dir"JCEC.raw.input.SE.txt"
echo "JCEC_RNASeq_SE_MATS="
echo $JCEC_RNASeq_SE_MATS

JC_RNASeq_A3SS_MATS_nohead=$dir"jcrnaseqa3ssmatsnohead.txt"
echo "JC_RNASeq_A3SS_MATS_nohead="
echo $JC_RNASeq_A3SS_MATS_nohead
JC_RNASeq_A5SS_MATS_nohead=$dir"jcrnaseqa5ssmatsnohead.txt"
echo "JC_RNASeq_A5SS_MATS_nohead="
echo $JC_RNASeq_A5SS_MATS_nohead
JC_RNASeq_MXE_MATS_nohead=$dir"jcrnaseqmxematsnohead.txt"
echo "JC_RNASeq_MXE_MATS_nohead="
echo $JC_RNASeq_MXE_MATS_nohead
JC_RNASeq_RI_MATS_nohead=$dir"jcrnaseqrimatsnohead.txt"
echo "JC_RNASeq_RI_MATS_nohead="
echo $JC_RNASeq_RI_MATS_nohead
JC_RNASeq_SE_MATS_nohead=$dir"jcrnaseqsematsnohead.txt"
echo "JC_RNASeq_SE_MATS_nohead="
echo $JC_RNASeq_SE_MATS_nohead

JCEC_RNASeq_A3SS_MATS_nohead=$dir"jcecrnaseqa3ssmatsnohead.txt"
echo "JCEC_RNASeq_A3SS_MATS_nohead="
echo $JCEC_RNASeq_A3SS_MATS_nohead
JCEC_RNASeq_A5SS_MATS_nohead=$dir"jcecrnaseqa5ssmatsnohead.txt"
echo "JCEC_RNASeq_A5SS_MATS_nohead="
echo $JCEC_RNASeq_A5SS_MATS_nohead
JCEC_RNASeq_MXE_MATS_nohead=$dir"jcecrnaseqmxematsnohead.txt"
echo "JCEC_RNASeq_MXE_MATS_nohead="
echo $JCEC_RNASeq_MXE_MATS_nohead
JCEC_RNASeq_RI_MATS_nohead=$dir"jcecrnaseqrimatsnohead.txt"
echo "JCEC_RNASeq_RI_MATS_nohead="
echo $JCEC_RNASeq_RI_MATS_nohead
JCEC_RNASeq_SE_MATS_nohead=$dir"jcecrnaseqsematsnohead.txt"
echo "JCEC_RNASeq_SE_MATS_nohead="
echo $JCEC_RNASeq_SE_MATS_nohead
#------------------------------------------------------------------------------------------------------
#
# Remove first line (column header) of each of the files
#
#------------------------------------------------------------------------------------------------------
sed '1d' $JC_RNASeq_A3SS_MATS > $JC_RNASeq_A3SS_MATS_nohead
sed '1d' $JC_RNASeq_A5SS_MATS > $JC_RNASeq_A5SS_MATS_nohead
sed '1d' $JC_RNASeq_MXE_MATS  > $JC_RNASeq_MXE_MATS_nohead
sed '1d' $JC_RNASeq_RI_MATS   > $JC_RNASeq_RI_MATS_nohead
sed '1d' $JC_RNASeq_SE_MATS   > $JC_RNASeq_SE_MATS_nohead

sed '1d' $JCEC_RNASeq_A3SS_MATS > $JCEC_RNASeq_A3SS_MATS_nohead
sed '1d' $JCEC_RNASeq_A5SS_MATS > $JCEC_RNASeq_A5SS_MATS_nohead
sed '1d' $JCEC_RNASeq_MXE_MATS  > $JCEC_RNASeq_MXE_MATS_nohead
sed '1d' $JCEC_RNASeq_RI_MATS   > $JCEC_RNASeq_RI_MATS_nohead
sed '1d' $JCEC_RNASeq_SE_MATS   > $JCEC_RNASeq_SE_MATS_nohead
#------------------------------------------------------------------------------------------------------
#
# Calculate inclusion level -- only difference between inclusion level calculations 
#                              for each of the different AS events is the
#                              IncFormLen and SkipFormLen
#                              We are accepting here those as calulated by rMATS
#
#------------------------------------------------------------------------------------------------------
#
# A3SS
#
#------------------------------------------------------------------------------------------------------
awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_A3SS_MATS_nohead   > $a3ss_jc_incend_sample1
awk -F "\t" '{print $1}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_idend_sample1
awk -F "\t" '{print $2}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_ijcend_sample1
awk -F "\t" '{print $3}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_sjcend_sample1
awk -F "\t" '{print $6}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_inclenend_sample1
awk -F "\t" '{print $7}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_skiplenend_sample1

awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_incend_sample1
awk -F "\t" '{print $1}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_idend_sample1
awk -F "\t" '{print $2}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_ijcend_sample1
awk -F "\t" '{print $3}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_sjcend_sample1
awk -F "\t" '{print $6}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_inclenend_sample1
awk -F "\t" '{print $7}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_skiplenend_sample1

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_A3SS_MATS_nohead   > $a3ss_jc_incend_sample2
awk -F "\t" '{print $1}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_idend_sample2
awk -F "\t" '{print $4}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_ijcend_sample2
awk -F "\t" '{print $5}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_sjcend_sample2
awk -F "\t" '{print $6}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_inclenend_sample2
awk -F "\t" '{print $7}' $JC_RNASeq_A3SS_MATS_nohead > $a3ss_jc_skiplenend_sample2

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_incend_sample2
awk -F "\t" '{print $1}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_idend_sample2
awk -F "\t" '{print $4}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_ijcend_sample2
awk -F "\t" '{print $5}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_sjcend_sample2
awk -F "\t" '{print $6}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_inclenend_sample2
awk -F "\t" '{print $7}' $JCEC_RNASeq_A3SS_MATS_nohead > $a3ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# A5SS
#
#------------------------------------------------------------------------------------------------------
awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_A5SS_MATS_nohead   > $a5ss_jc_incend_sample1
awk -F "\t" '{print $1}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_idend_sample1
awk -F "\t" '{print $2}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_ijcend_sample1
awk -F "\t" '{print $3}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_sjcend_sample1
awk -F "\t" '{print $6}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_inclenend_sample1
awk -F "\t" '{print $7}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_skiplenend_sample1

awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_incend_sample1
awk -F "\t" '{print $1}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_idend_sample1
awk -F "\t" '{print $2}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_ijcend_sample1
awk -F "\t" '{print $3}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_sjcend_sample1
awk -F "\t" '{print $6}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_inclenend_sample1
awk -F "\t" '{print $7}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_skiplenend_sample1

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_A5SS_MATS_nohead   > $a5ss_jc_incend_sample2
awk -F "\t" '{print $1}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_idend_sample2
awk -F "\t" '{print $4}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_ijcend_sample2
awk -F "\t" '{print $5}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_sjcend_sample2
awk -F "\t" '{print $6}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_inclenend_sample2
awk -F "\t" '{print $7}' $JC_RNASeq_A5SS_MATS_nohead > $a5ss_jc_skiplenend_sample2

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_incend_sample2
awk -F "\t" '{print $1}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_idend_sample2
awk -F "\t" '{print $4}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_ijcend_sample2
awk -F "\t" '{print $5}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_sjcend_sample2
awk -F "\t" '{print $6}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_inclenend_sample2
awk -F "\t" '{print $7}' $JCEC_RNASeq_A5SS_MATS_nohead > $a5ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# MXE
#
#------------------------------------------------------------------------------------------------------
awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_MXE_MATS_nohead   > $mxe_jc_incend_sample1
awk -F "\t" '{print $1}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_idend_sample1
awk -F "\t" '{print $2}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_ijcend_sample1
awk -F "\t" '{print $3}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_sjcend_sample1
awk -F "\t" '{print $6}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_inclenend_sample1
awk -F "\t" '{print $7}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_skiplenend_sample1

awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_incend_sample1
awk -F "\t" '{print $1}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_idend_sample1
awk -F "\t" '{print $2}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_ijcend_sample1
awk -F "\t" '{print $3}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_sjcend_sample1
awk -F "\t" '{print $6}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_inclenend_sample1
awk -F "\t" '{print $7}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_skiplenend_sample1

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_MXE_MATS_nohead   > $mxe_jc_incend_sample2
awk -F "\t" '{print $1}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_idend_sample2
awk -F "\t" '{print $4}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_ijcend_sample2
awk -F "\t" '{print $5}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_sjcend_sample2
awk -F "\t" '{print $6}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_inclenend_sample2
awk -F "\t" '{print $7}' $JC_RNASeq_MXE_MATS_nohead > $mxe_jc_skiplenend_sample2

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_incend_sample2
awk -F "\t" '{print $1}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_idend_sample2
awk -F "\t" '{print $4}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_ijcend_sample2
awk -F "\t" '{print $5}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_sjcend_sample2
awk -F "\t" '{print $6}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_inclenend_sample2
awk -F "\t" '{print $7}' $JCEC_RNASeq_MXE_MATS_nohead > $mxe_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# RI
#
#------------------------------------------------------------------------------------------------------
awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_RI_MATS_nohead   > $ri_jc_incend_sample1
awk -F "\t" '{print $1}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_idend_sample1
awk -F "\t" '{print $2}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_ijcend_sample1
awk -F "\t" '{print $3}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_sjcend_sample1
awk -F "\t" '{print $6}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_inclenend_sample1
awk -F "\t" '{print $7}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_skiplenend_sample1

awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_incend_sample1
awk -F "\t" '{print $1}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_idend_sample1
awk -F "\t" '{print $2}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_ijcend_sample1
awk -F "\t" '{print $3}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_sjcend_sample1
awk -F "\t" '{print $6}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_inclenend_sample1
awk -F "\t" '{print $7}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_skiplenend_sample1

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_RI_MATS_nohead   > $ri_jc_incend_sample2
awk -F "\t" '{print $1}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_idend_sample2
awk -F "\t" '{print $4}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_ijcend_sample2
awk -F "\t" '{print $5}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_sjcend_sample2
awk -F "\t" '{print $6}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_inclenend_sample2
awk -F "\t" '{print $7}' $JC_RNASeq_RI_MATS_nohead > $ri_jc_skiplenend_sample2

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_incend_sample2
awk -F "\t" '{print $1}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_idend_sample2
awk -F "\t" '{print $4}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_ijcend_sample2
awk -F "\t" '{print $5}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_sjcend_sample2
awk -F "\t" '{print $6}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_inclenend_sample2
awk -F "\t" '{print $7}' $JCEC_RNASeq_RI_MATS_nohead > $ri_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# SE
#
#------------------------------------------------------------------------------------------------------
awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_SE_MATS_nohead   > $se_jc_incend_sample1
awk -F "\t" '{print $1}' $JC_RNASeq_SE_MATS_nohead > $se_jc_idend_sample1
awk -F "\t" '{print $2}' $JC_RNASeq_SE_MATS_nohead > $se_jc_ijcend_sample1
awk -F "\t" '{print $3}' $JC_RNASeq_SE_MATS_nohead > $se_jc_sjcend_sample1
awk -F "\t" '{print $6}' $JC_RNASeq_SE_MATS_nohead > $se_jc_inclenend_sample1
awk -F "\t" '{print $7}' $JC_RNASeq_SE_MATS_nohead > $se_jc_skiplenend_sample1

awk -F "\t" '{num=$2/$6; den=($2/$6) + ($3/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_incend_sample1
awk -F "\t" '{print $1}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_idend_sample1
awk -F "\t" '{print $2}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_ijcend_sample1
awk -F "\t" '{print $3}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_sjcend_sample1
awk -F "\t" '{print $6}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_inclenend_sample1
awk -F "\t" '{print $7}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_skiplenend_sample1

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JC_RNASeq_SE_MATS_nohead   > $se_jc_incend_sample2
awk -F "\t" '{print $1}' $JC_RNASeq_SE_MATS_nohead > $se_jc_idend_sample2
awk -F "\t" '{print $4}' $JC_RNASeq_SE_MATS_nohead > $se_jc_ijcend_sample2
awk -F "\t" '{print $5}' $JC_RNASeq_SE_MATS_nohead > $se_jc_sjcend_sample2
awk -F "\t" '{print $6}' $JC_RNASeq_SE_MATS_nohead > $se_jc_inclenend_sample2
awk -F "\t" '{print $7}' $JC_RNASeq_SE_MATS_nohead > $se_jc_skiplenend_sample2

awk -F "\t" '{num=$4/$6; den=($4/$6) + ($5/$7); if (den > 0) {inc=num/den;} else {inc=0;} printf "%.6s\n", inc;}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_incend_sample2
awk -F "\t" '{print $1}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_idend_sample2
awk -F "\t" '{print $4}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_ijcend_sample2
awk -F "\t" '{print $5}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_sjcend_sample2
awk -F "\t" '{print $6}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_inclenend_sample2
awk -F "\t" '{print $7}' $JCEC_RNASeq_SE_MATS_nohead > $se_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# now prepend the sample name to the top of all these columns
#
#------------------------------------------------------------------------------------------------------
#
# A3SS
#
#------------------------------------------------------------------------------------------------------
sed "1 i\\ID"       $a3ss_jc_idend_sample1        > tmp; mv tmp $a3ss_jc_idend_sample1
sed "1 i\\$sample1" $a3ss_jc_incend_sample1       > tmp; mv tmp $a3ss_jc_incend_sample1
sed "1 i\\$sample1" $a3ss_jc_ijcend_sample1       > tmp; mv tmp $a3ss_jc_ijcend_sample1
sed "1 i\\$sample1" $a3ss_jc_sjcend_sample1       > tmp; mv tmp $a3ss_jc_sjcend_sample1
sed "1 i\\$sample1" $a3ss_jc_inclenend_sample1    > tmp; mv tmp $a3ss_jc_inclenend_sample1
sed "1 i\\$sample1" $a3ss_jc_skiplenend_sample1   > tmp; mv tmp $a3ss_jc_skiplenend_sample1

sed "1 i\\ID"       $a3ss_jcec_idend_sample1      > tmp; mv tmp $a3ss_jcec_idend_sample1
sed "1 i\\$sample1" $a3ss_jcec_incend_sample1     > tmp; mv tmp $a3ss_jcec_incend_sample1
sed "1 i\\$sample1" $a3ss_jcec_ijcend_sample1     > tmp; mv tmp $a3ss_jcec_ijcend_sample1
sed "1 i\\$sample1" $a3ss_jcec_sjcend_sample1     > tmp; mv tmp $a3ss_jcec_sjcend_sample1
sed "1 i\\$sample1" $a3ss_jcec_inclenend_sample1  > tmp; mv tmp $a3ss_jcec_inclenend_sample1
sed "1 i\\$sample1" $a3ss_jcec_skiplenend_sample1 > tmp; mv tmp $a3ss_jcec_skiplenend_sample1

sed "1 i\\ID"       $a3ss_jc_idend_sample2        > tmp; mv tmp $a3ss_jc_idend_sample2
sed "1 i\\$sample2" $a3ss_jc_incend_sample2       > tmp; mv tmp $a3ss_jc_incend_sample2
sed "1 i\\$sample2" $a3ss_jc_ijcend_sample2       > tmp; mv tmp $a3ss_jc_ijcend_sample2
sed "1 i\\$sample2" $a3ss_jc_sjcend_sample2       > tmp; mv tmp $a3ss_jc_sjcend_sample2
sed "1 i\\$sample2" $a3ss_jc_inclenend_sample2    > tmp; mv tmp $a3ss_jc_inclenend_sample2
sed "1 i\\$sample2" $a3ss_jc_skiplenend_sample2   > tmp; mv tmp $a3ss_jc_skiplenend_sample2

sed "1 i\\ID"       $a3ss_jcec_idend_sample2      > tmp; mv tmp $a3ss_jcec_idend_sample2
sed "1 i\\$sample2" $a3ss_jcec_incend_sample2     > tmp; mv tmp $a3ss_jcec_incend_sample2
sed "1 i\\$sample2" $a3ss_jcec_ijcend_sample2     > tmp; mv tmp $a3ss_jcec_ijcend_sample2
sed "1 i\\$sample2" $a3ss_jcec_sjcend_sample2     > tmp; mv tmp $a3ss_jcec_sjcend_sample2
sed "1 i\\$sample2" $a3ss_jcec_inclenend_sample2  > tmp; mv tmp $a3ss_jcec_inclenend_sample2
sed "1 i\\$sample2" $a3ss_jcec_skiplenend_sample2 > tmp; mv tmp $a3ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# A5SS
#
#------------------------------------------------------------------------------------------------------
sed "1 i\\ID"       $a5ss_jc_idend_sample1        > tmp; mv tmp $a5ss_jc_idend_sample1
sed "1 i\\$sample1" $a5ss_jc_incend_sample1       > tmp; mv tmp $a5ss_jc_incend_sample1
sed "1 i\\$sample1" $a5ss_jc_ijcend_sample1       > tmp; mv tmp $a5ss_jc_ijcend_sample1
sed "1 i\\$sample1" $a5ss_jc_sjcend_sample1       > tmp; mv tmp $a5ss_jc_sjcend_sample1
sed "1 i\\$sample1" $a5ss_jc_inclenend_sample1    > tmp; mv tmp $a5ss_jc_inclenend_sample1
sed "1 i\\$sample1" $a5ss_jc_skiplenend_sample1   > tmp; mv tmp $a5ss_jc_skiplenend_sample1

sed "1 i\\ID"       $a5ss_jcec_idend_sample1      > tmp; mv tmp $a5ss_jcec_idend_sample1
sed "1 i\\$sample1" $a5ss_jcec_incend_sample1     > tmp; mv tmp $a5ss_jcec_incend_sample1
sed "1 i\\$sample1" $a5ss_jcec_ijcend_sample1     > tmp; mv tmp $a5ss_jcec_ijcend_sample1
sed "1 i\\$sample1" $a5ss_jcec_sjcend_sample1     > tmp; mv tmp $a5ss_jcec_sjcend_sample1
sed "1 i\\$sample1" $a5ss_jcec_inclenend_sample1  > tmp; mv tmp $a5ss_jcec_inclenend_sample1
sed "1 i\\$sample1" $a5ss_jcec_skiplenend_sample1 > tmp; mv tmp $a5ss_jcec_skiplenend_sample1

sed "1 i\\ID"       $a5ss_jc_idend_sample2        > tmp; mv tmp $a5ss_jc_idend_sample2
sed "1 i\\$sample2" $a5ss_jc_incend_sample2       > tmp; mv tmp $a5ss_jc_incend_sample2
sed "1 i\\$sample2" $a5ss_jc_ijcend_sample2       > tmp; mv tmp $a5ss_jc_ijcend_sample2
sed "1 i\\$sample2" $a5ss_jc_sjcend_sample2       > tmp; mv tmp $a5ss_jc_sjcend_sample2
sed "1 i\\$sample2" $a5ss_jc_inclenend_sample2    > tmp; mv tmp $a5ss_jc_inclenend_sample2
sed "1 i\\$sample2" $a5ss_jc_skiplenend_sample2   > tmp; mv tmp $a5ss_jc_skiplenend_sample2

sed "1 i\\ID"       $a5ss_jcec_idend_sample2      > tmp; mv tmp $a5ss_jcec_idend_sample2
sed "1 i\\$sample2" $a5ss_jcec_incend_sample2     > tmp; mv tmp $a5ss_jcec_incend_sample2
sed "1 i\\$sample2" $a5ss_jcec_ijcend_sample2     > tmp; mv tmp $a5ss_jcec_ijcend_sample2
sed "1 i\\$sample2" $a5ss_jcec_sjcend_sample2     > tmp; mv tmp $a5ss_jcec_sjcend_sample2
sed "1 i\\$sample2" $a5ss_jcec_inclenend_sample2  > tmp; mv tmp $a5ss_jcec_inclenend_sample2
sed "1 i\\$sample2" $a5ss_jcec_skiplenend_sample2 > tmp; mv tmp $a5ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# MXE
#
#------------------------------------------------------------------------------------------------------
sed "1 i\\ID"       $mxe_jc_idend_sample1        > tmp; mv tmp $mxe_jc_idend_sample1
sed "1 i\\$sample1" $mxe_jc_incend_sample1       > tmp; mv tmp $mxe_jc_incend_sample1
sed "1 i\\$sample1" $mxe_jc_ijcend_sample1       > tmp; mv tmp $mxe_jc_ijcend_sample1
sed "1 i\\$sample1" $mxe_jc_sjcend_sample1       > tmp; mv tmp $mxe_jc_sjcend_sample1
sed "1 i\\$sample1" $mxe_jc_inclenend_sample1    > tmp; mv tmp $mxe_jc_inclenend_sample1
sed "1 i\\$sample1" $mxe_jc_skiplenend_sample1   > tmp; mv tmp $mxe_jc_skiplenend_sample1

sed "1 i\\ID"       $mxe_jcec_idend_sample1      > tmp; mv tmp $mxe_jcec_idend_sample1
sed "1 i\\$sample1" $mxe_jcec_incend_sample1     > tmp; mv tmp $mxe_jcec_incend_sample1
sed "1 i\\$sample1" $mxe_jcec_ijcend_sample1     > tmp; mv tmp $mxe_jcec_ijcend_sample1
sed "1 i\\$sample1" $mxe_jcec_sjcend_sample1     > tmp; mv tmp $mxe_jcec_sjcend_sample1
sed "1 i\\$sample1" $mxe_jcec_inclenend_sample1  > tmp; mv tmp $mxe_jcec_inclenend_sample1
sed "1 i\\$sample1" $mxe_jcec_skiplenend_sample1 > tmp; mv tmp $mxe_jcec_skiplenend_sample1

sed "1 i\\ID"       $mxe_jc_idend_sample2        > tmp; mv tmp $mxe_jc_idend_sample2
sed "1 i\\$sample2" $mxe_jc_incend_sample2       > tmp; mv tmp $mxe_jc_incend_sample2
sed "1 i\\$sample2" $mxe_jc_ijcend_sample2       > tmp; mv tmp $mxe_jc_ijcend_sample2
sed "1 i\\$sample2" $mxe_jc_sjcend_sample2       > tmp; mv tmp $mxe_jc_sjcend_sample2
sed "1 i\\$sample2" $mxe_jc_inclenend_sample2    > tmp; mv tmp $mxe_jc_inclenend_sample2
sed "1 i\\$sample2" $mxe_jc_skiplenend_sample2   > tmp; mv tmp $mxe_jc_skiplenend_sample2

sed "1 i\\ID"       $mxe_jcec_idend_sample2      > tmp; mv tmp $mxe_jcec_idend_sample2
sed "1 i\\$sample2" $mxe_jcec_incend_sample2     > tmp; mv tmp $mxe_jcec_incend_sample2
sed "1 i\\$sample2" $mxe_jcec_ijcend_sample2     > tmp; mv tmp $mxe_jcec_ijcend_sample2
sed "1 i\\$sample2" $mxe_jcec_sjcend_sample2     > tmp; mv tmp $mxe_jcec_sjcend_sample2
sed "1 i\\$sample2" $mxe_jcec_inclenend_sample2  > tmp; mv tmp $mxe_jcec_inclenend_sample2
sed "1 i\\$sample2" $mxe_jcec_skiplenend_sample2 > tmp; mv tmp $mxe_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# RI
#
#------------------------------------------------------------------------------------------------------
sed "1 i\\ID"       $ri_jc_idend_sample1        > tmp; mv tmp $ri_jc_idend_sample1
sed "1 i\\$sample1" $ri_jc_incend_sample1       > tmp; mv tmp $ri_jc_incend_sample1
sed "1 i\\$sample1" $ri_jc_ijcend_sample1       > tmp; mv tmp $ri_jc_ijcend_sample1
sed "1 i\\$sample1" $ri_jc_sjcend_sample1       > tmp; mv tmp $ri_jc_sjcend_sample1
sed "1 i\\$sample1" $ri_jc_inclenend_sample1    > tmp; mv tmp $ri_jc_inclenend_sample1
sed "1 i\\$sample1" $ri_jc_skiplenend_sample1   > tmp; mv tmp $ri_jc_skiplenend_sample1

sed "1 i\\ID"       $ri_jcec_idend_sample1      > tmp; mv tmp $ri_jcec_idend_sample1
sed "1 i\\$sample1" $ri_jcec_incend_sample1     > tmp; mv tmp $ri_jcec_incend_sample1
sed "1 i\\$sample1" $ri_jcec_ijcend_sample1     > tmp; mv tmp $ri_jcec_ijcend_sample1
sed "1 i\\$sample1" $ri_jcec_sjcend_sample1     > tmp; mv tmp $ri_jcec_sjcend_sample1
sed "1 i\\$sample1" $ri_jcec_inclenend_sample1  > tmp; mv tmp $ri_jcec_inclenend_sample1
sed "1 i\\$sample1" $ri_jcec_skiplenend_sample1 > tmp; mv tmp $ri_jcec_skiplenend_sample1

sed "1 i\\ID"       $ri_jc_idend_sample2        > tmp; mv tmp $ri_jc_idend_sample2
sed "1 i\\$sample2" $ri_jc_incend_sample2       > tmp; mv tmp $ri_jc_incend_sample2
sed "1 i\\$sample2" $ri_jc_ijcend_sample2       > tmp; mv tmp $ri_jc_ijcend_sample2
sed "1 i\\$sample2" $ri_jc_sjcend_sample2       > tmp; mv tmp $ri_jc_sjcend_sample2
sed "1 i\\$sample2" $ri_jc_inclenend_sample2    > tmp; mv tmp $ri_jc_inclenend_sample2
sed "1 i\\$sample2" $ri_jc_skiplenend_sample2   > tmp; mv tmp $ri_jc_skiplenend_sample2

sed "1 i\\ID"       $ri_jcec_idend_sample2      > tmp; mv tmp $ri_jcec_idend_sample2
sed "1 i\\$sample2" $ri_jcec_incend_sample2     > tmp; mv tmp $ri_jcec_incend_sample2
sed "1 i\\$sample2" $ri_jcec_ijcend_sample2     > tmp; mv tmp $ri_jcec_ijcend_sample2
sed "1 i\\$sample2" $ri_jcec_sjcend_sample2     > tmp; mv tmp $ri_jcec_sjcend_sample2
sed "1 i\\$sample2" $ri_jcec_inclenend_sample2  > tmp; mv tmp $ri_jcec_inclenend_sample2
sed "1 i\\$sample2" $ri_jcec_skiplenend_sample2 > tmp; mv tmp $ri_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# SE
#
#------------------------------------------------------------------------------------------------------
sed "1 i\\ID"       $se_jc_idend_sample1        > tmp; mv tmp $se_jc_idend_sample1
sed "1 i\\$sample1" $se_jc_incend_sample1       > tmp; mv tmp $se_jc_incend_sample1
sed "1 i\\$sample1" $se_jc_ijcend_sample1       > tmp; mv tmp $se_jc_ijcend_sample1
sed "1 i\\$sample1" $se_jc_sjcend_sample1       > tmp; mv tmp $se_jc_sjcend_sample1
sed "1 i\\$sample1" $se_jc_inclenend_sample1    > tmp; mv tmp $se_jc_inclenend_sample1
sed "1 i\\$sample1" $se_jc_skiplenend_sample1   > tmp; mv tmp $se_jc_skiplenend_sample1

sed "1 i\\ID"       $se_jcec_idend_sample1      > tmp; mv tmp $se_jcec_idend_sample1
sed "1 i\\$sample1" $se_jcec_incend_sample1     > tmp; mv tmp $se_jcec_incend_sample1
sed "1 i\\$sample1" $se_jcec_ijcend_sample1     > tmp; mv tmp $se_jcec_ijcend_sample1
sed "1 i\\$sample1" $se_jcec_sjcend_sample1     > tmp; mv tmp $se_jcec_sjcend_sample1
sed "1 i\\$sample1" $se_jcec_inclenend_sample1  > tmp; mv tmp $se_jcec_inclenend_sample1
sed "1 i\\$sample1" $se_jcec_skiplenend_sample1 > tmp; mv tmp $se_jcec_skiplenend_sample1

sed "1 i\\ID"       $se_jc_idend_sample2        > tmp; mv tmp $se_jc_idend_sample2
sed "1 i\\$sample2" $se_jc_incend_sample2       > tmp; mv tmp $se_jc_incend_sample2
sed "1 i\\$sample2" $se_jc_ijcend_sample2       > tmp; mv tmp $se_jc_ijcend_sample2
sed "1 i\\$sample2" $se_jc_sjcend_sample2       > tmp; mv tmp $se_jc_sjcend_sample2
sed "1 i\\$sample2" $se_jc_inclenend_sample2    > tmp; mv tmp $se_jc_inclenend_sample2
sed "1 i\\$sample2" $se_jc_skiplenend_sample2   > tmp; mv tmp $se_jc_skiplenend_sample2

sed "1 i\\ID"       $se_jcec_idend_sample2      > tmp; mv tmp $se_jcec_idend_sample2
sed "1 i\\$sample2" $se_jcec_incend_sample2     > tmp; mv tmp $se_jcec_incend_sample2
sed "1 i\\$sample2" $se_jcec_ijcend_sample2     > tmp; mv tmp $se_jcec_ijcend_sample2
sed "1 i\\$sample2" $se_jcec_sjcend_sample2     > tmp; mv tmp $se_jcec_sjcend_sample2
sed "1 i\\$sample2" $se_jcec_inclenend_sample2  > tmp; mv tmp $se_jcec_inclenend_sample2
sed "1 i\\$sample2" $se_jcec_skiplenend_sample2 > tmp; mv tmp $se_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# make sample specific matrices with id with $sample1 and sample2
#------------------------------------------------------------------------------------------------------
#
# A3SS 
#
#------------------------------------------------------------------------------------------------------
paste -d ' ' $a3ss_jc_idend_sample1 $a3ss_jc_incend_sample1         > tmp; mv tmp $a3ss_jc_incend_sample1
paste -d ' ' $a3ss_jc_idend_sample1 $a3ss_jc_ijcend_sample1         > tmp; mv tmp $a3ss_jc_ijcend_sample1
paste -d ' ' $a3ss_jc_idend_sample1 $a3ss_jc_sjcend_sample1         > tmp; mv tmp $a3ss_jc_sjcend_sample1
paste -d ' ' $a3ss_jc_idend_sample1 $a3ss_jc_inclenend_sample1      > tmp; mv tmp $a3ss_jc_inclenend_sample1
paste -d ' ' $a3ss_jc_idend_sample1 $a3ss_jc_skiplenend_sample1     > tmp; mv tmp $a3ss_jc_skiplenend_sample1

paste -d ' ' $a3ss_jcec_idend_sample1 $a3ss_jcec_incend_sample1     > tmp; mv tmp $a3ss_jcec_incend_sample1
paste -d ' ' $a3ss_jcec_idend_sample1 $a3ss_jcec_ijcend_sample1     > tmp; mv tmp $a3ss_jcec_ijcend_sample1
paste -d ' ' $a3ss_jcec_idend_sample1 $a3ss_jcec_sjcend_sample1     > tmp; mv tmp $a3ss_jcec_sjcend_sample1
paste -d ' ' $a3ss_jcec_idend_sample1 $a3ss_jcec_inclenend_sample1  > tmp; mv tmp $a3ss_jcec_inclenend_sample1
paste -d ' ' $a3ss_jcec_idend_sample1 $a3ss_jcec_skiplenend_sample1 > tmp; mv tmp $a3ss_jcec_skiplenend_sample1

paste -d ' ' $a3ss_jc_idend_sample2 $a3ss_jc_incend_sample2         > tmp; mv tmp $a3ss_jc_incend_sample2
paste -d ' ' $a3ss_jc_idend_sample2 $a3ss_jc_ijcend_sample2         > tmp; mv tmp $a3ss_jc_ijcend_sample2
paste -d ' ' $a3ss_jc_idend_sample2 $a3ss_jc_sjcend_sample2         > tmp; mv tmp $a3ss_jc_sjcend_sample2
paste -d ' ' $a3ss_jc_idend_sample2 $a3ss_jc_inclenend_sample2      > tmp; mv tmp $a3ss_jc_inclenend_sample2
paste -d ' ' $a3ss_jc_idend_sample2 $a3ss_jc_skiplenend_sample2     > tmp; mv tmp $a3ss_jc_skiplenend_sample2

paste -d ' ' $a3ss_jcec_idend_sample2 $a3ss_jcec_incend_sample2     > tmp; mv tmp $a3ss_jcec_incend_sample2
paste -d ' ' $a3ss_jcec_idend_sample2 $a3ss_jcec_ijcend_sample2     > tmp; mv tmp $a3ss_jcec_ijcend_sample2
paste -d ' ' $a3ss_jcec_idend_sample2 $a3ss_jcec_sjcend_sample2     > tmp; mv tmp $a3ss_jcec_sjcend_sample2
paste -d ' ' $a3ss_jcec_idend_sample2 $a3ss_jcec_inclenend_sample2  > tmp; mv tmp $a3ss_jcec_inclenend_sample2
paste -d ' ' $a3ss_jcec_idend_sample2 $a3ss_jcec_skiplenend_sample2 > tmp; mv tmp $a3ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# A5SS 
#
#------------------------------------------------------------------------------------------------------
paste -d ' ' $a5ss_jc_idend_sample1 $a5ss_jc_incend_sample1         > tmp; mv tmp $a5ss_jc_incend_sample1
paste -d ' ' $a5ss_jc_idend_sample1 $a5ss_jc_ijcend_sample1         > tmp; mv tmp $a5ss_jc_ijcend_sample1
paste -d ' ' $a5ss_jc_idend_sample1 $a5ss_jc_sjcend_sample1         > tmp; mv tmp $a5ss_jc_sjcend_sample1
paste -d ' ' $a5ss_jc_idend_sample1 $a5ss_jc_inclenend_sample1      > tmp; mv tmp $a5ss_jc_inclenend_sample1
paste -d ' ' $a5ss_jc_idend_sample1 $a5ss_jc_skiplenend_sample1     > tmp; mv tmp $a5ss_jc_skiplenend_sample1

paste -d ' ' $a5ss_jcec_idend_sample1 $a5ss_jcec_incend_sample1     > tmp; mv tmp $a5ss_jcec_incend_sample1
paste -d ' ' $a5ss_jcec_idend_sample1 $a5ss_jcec_ijcend_sample1     > tmp; mv tmp $a5ss_jcec_ijcend_sample1
paste -d ' ' $a5ss_jcec_idend_sample1 $a5ss_jcec_sjcend_sample1     > tmp; mv tmp $a5ss_jcec_sjcend_sample1
paste -d ' ' $a5ss_jcec_idend_sample1 $a5ss_jcec_inclenend_sample1  > tmp; mv tmp $a5ss_jcec_inclenend_sample1
paste -d ' ' $a5ss_jcec_idend_sample1 $a5ss_jcec_skiplenend_sample1 > tmp; mv tmp $a5ss_jcec_skiplenend_sample1

paste -d ' ' $a5ss_jc_idend_sample2 $a5ss_jc_incend_sample2         > tmp; mv tmp $a5ss_jc_incend_sample2
paste -d ' ' $a5ss_jc_idend_sample2 $a5ss_jc_ijcend_sample2         > tmp; mv tmp $a5ss_jc_ijcend_sample2
paste -d ' ' $a5ss_jc_idend_sample2 $a5ss_jc_sjcend_sample2         > tmp; mv tmp $a5ss_jc_sjcend_sample2
paste -d ' ' $a5ss_jc_idend_sample2 $a5ss_jc_inclenend_sample2      > tmp; mv tmp $a5ss_jc_inclenend_sample2
paste -d ' ' $a5ss_jc_idend_sample2 $a5ss_jc_skiplenend_sample2     > tmp; mv tmp $a5ss_jc_skiplenend_sample2

paste -d ' ' $a5ss_jcec_idend_sample2 $a5ss_jcec_incend_sample2     > tmp; mv tmp $a5ss_jcec_incend_sample2
paste -d ' ' $a5ss_jcec_idend_sample2 $a5ss_jcec_ijcend_sample2     > tmp; mv tmp $a5ss_jcec_ijcend_sample2
paste -d ' ' $a5ss_jcec_idend_sample2 $a5ss_jcec_sjcend_sample2     > tmp; mv tmp $a5ss_jcec_sjcend_sample2
paste -d ' ' $a5ss_jcec_idend_sample2 $a5ss_jcec_inclenend_sample2  > tmp; mv tmp $a5ss_jcec_inclenend_sample2
paste -d ' ' $a5ss_jcec_idend_sample2 $a5ss_jcec_skiplenend_sample2 > tmp; mv tmp $a5ss_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# MXE 
#
#------------------------------------------------------------------------------------------------------
paste -d ' ' $mxe_jc_idend_sample1 $mxe_jc_incend_sample1         > tmp; mv tmp $mxe_jc_incend_sample1
paste -d ' ' $mxe_jc_idend_sample1 $mxe_jc_ijcend_sample1         > tmp; mv tmp $mxe_jc_ijcend_sample1
paste -d ' ' $mxe_jc_idend_sample1 $mxe_jc_sjcend_sample1         > tmp; mv tmp $mxe_jc_sjcend_sample1
paste -d ' ' $mxe_jc_idend_sample1 $mxe_jc_inclenend_sample1      > tmp; mv tmp $mxe_jc_inclenend_sample1
paste -d ' ' $mxe_jc_idend_sample1 $mxe_jc_skiplenend_sample1     > tmp; mv tmp $mxe_jc_skiplenend_sample1

paste -d ' ' $mxe_jcec_idend_sample1 $mxe_jcec_incend_sample1     > tmp; mv tmp $mxe_jcec_incend_sample1
paste -d ' ' $mxe_jcec_idend_sample1 $mxe_jcec_ijcend_sample1     > tmp; mv tmp $mxe_jcec_ijcend_sample1
paste -d ' ' $mxe_jcec_idend_sample1 $mxe_jcec_sjcend_sample1     > tmp; mv tmp $mxe_jcec_sjcend_sample1
paste -d ' ' $mxe_jcec_idend_sample1 $mxe_jcec_inclenend_sample1  > tmp; mv tmp $mxe_jcec_inclenend_sample1
paste -d ' ' $mxe_jcec_idend_sample1 $mxe_jcec_skiplenend_sample1 > tmp; mv tmp $mxe_jcec_skiplenend_sample1

paste -d ' ' $mxe_jc_idend_sample2 $mxe_jc_incend_sample2         > tmp; mv tmp $mxe_jc_incend_sample2
paste -d ' ' $mxe_jc_idend_sample2 $mxe_jc_ijcend_sample2         > tmp; mv tmp $mxe_jc_ijcend_sample2
paste -d ' ' $mxe_jc_idend_sample2 $mxe_jc_sjcend_sample2         > tmp; mv tmp $mxe_jc_sjcend_sample2
paste -d ' ' $mxe_jc_idend_sample2 $mxe_jc_inclenend_sample2      > tmp; mv tmp $mxe_jc_inclenend_sample2
paste -d ' ' $mxe_jc_idend_sample2 $mxe_jc_skiplenend_sample2     > tmp; mv tmp $mxe_jc_skiplenend_sample2

paste -d ' ' $mxe_jcec_idend_sample2 $mxe_jcec_incend_sample2     > tmp; mv tmp $mxe_jcec_incend_sample2
paste -d ' ' $mxe_jcec_idend_sample2 $mxe_jcec_ijcend_sample2     > tmp; mv tmp $mxe_jcec_ijcend_sample2
paste -d ' ' $mxe_jcec_idend_sample2 $mxe_jcec_sjcend_sample2     > tmp; mv tmp $mxe_jcec_sjcend_sample2
paste -d ' ' $mxe_jcec_idend_sample2 $mxe_jcec_inclenend_sample2  > tmp; mv tmp $mxe_jcec_inclenend_sample2
paste -d ' ' $mxe_jcec_idend_sample2 $mxe_jcec_skiplenend_sample2 > tmp; mv tmp $mxe_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# RI 
#
#------------------------------------------------------------------------------------------------------
paste -d ' ' $ri_jc_idend_sample1 $ri_jc_incend_sample1         > tmp; mv tmp $ri_jc_incend_sample1
paste -d ' ' $ri_jc_idend_sample1 $ri_jc_ijcend_sample1         > tmp; mv tmp $ri_jc_ijcend_sample1
paste -d ' ' $ri_jc_idend_sample1 $ri_jc_sjcend_sample1         > tmp; mv tmp $ri_jc_sjcend_sample1
paste -d ' ' $ri_jc_idend_sample1 $ri_jc_inclenend_sample1      > tmp; mv tmp $ri_jc_inclenend_sample1
paste -d ' ' $ri_jc_idend_sample1 $ri_jc_skiplenend_sample1     > tmp; mv tmp $ri_jc_skiplenend_sample1

paste -d ' ' $ri_jcec_idend_sample1 $ri_jcec_incend_sample1     > tmp; mv tmp $ri_jcec_incend_sample1
paste -d ' ' $ri_jcec_idend_sample1 $ri_jcec_ijcend_sample1     > tmp; mv tmp $ri_jcec_ijcend_sample1
paste -d ' ' $ri_jcec_idend_sample1 $ri_jcec_sjcend_sample1     > tmp; mv tmp $ri_jcec_sjcend_sample1
paste -d ' ' $ri_jcec_idend_sample1 $ri_jcec_inclenend_sample1  > tmp; mv tmp $ri_jcec_inclenend_sample1
paste -d ' ' $ri_jcec_idend_sample1 $ri_jcec_skiplenend_sample1 > tmp; mv tmp $ri_jcec_skiplenend_sample1

paste -d ' ' $ri_jc_idend_sample2 $ri_jc_incend_sample2         > tmp; mv tmp $ri_jc_incend_sample2
paste -d ' ' $ri_jc_idend_sample2 $ri_jc_ijcend_sample2         > tmp; mv tmp $ri_jc_ijcend_sample2
paste -d ' ' $ri_jc_idend_sample2 $ri_jc_sjcend_sample2         > tmp; mv tmp $ri_jc_sjcend_sample2
paste -d ' ' $ri_jc_idend_sample2 $ri_jc_inclenend_sample2      > tmp; mv tmp $ri_jc_inclenend_sample2
paste -d ' ' $ri_jc_idend_sample2 $ri_jc_skiplenend_sample2     > tmp; mv tmp $ri_jc_skiplenend_sample2

paste -d ' ' $ri_jcec_idend_sample2 $ri_jcec_incend_sample2     > tmp; mv tmp $ri_jcec_incend_sample2
paste -d ' ' $ri_jcec_idend_sample2 $ri_jcec_ijcend_sample2     > tmp; mv tmp $ri_jcec_ijcend_sample2
paste -d ' ' $ri_jcec_idend_sample2 $ri_jcec_sjcend_sample2     > tmp; mv tmp $ri_jcec_sjcend_sample2
paste -d ' ' $ri_jcec_idend_sample2 $ri_jcec_inclenend_sample2  > tmp; mv tmp $ri_jcec_inclenend_sample2
paste -d ' ' $ri_jcec_idend_sample2 $ri_jcec_skiplenend_sample2 > tmp; mv tmp $ri_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# SE 
#
#------------------------------------------------------------------------------------------------------
paste -d ' ' $se_jc_idend_sample1 $se_jc_incend_sample1         > tmp; mv tmp $se_jc_incend_sample1
paste -d ' ' $se_jc_idend_sample1 $se_jc_ijcend_sample1         > tmp; mv tmp $se_jc_ijcend_sample1
paste -d ' ' $se_jc_idend_sample1 $se_jc_sjcend_sample1         > tmp; mv tmp $se_jc_sjcend_sample1
paste -d ' ' $se_jc_idend_sample1 $se_jc_inclenend_sample1      > tmp; mv tmp $se_jc_inclenend_sample1
paste -d ' ' $se_jc_idend_sample1 $se_jc_skiplenend_sample1     > tmp; mv tmp $se_jc_skiplenend_sample1

paste -d ' ' $se_jcec_idend_sample1 $se_jcec_incend_sample1     > tmp; mv tmp $se_jcec_incend_sample1
paste -d ' ' $se_jcec_idend_sample1 $se_jcec_ijcend_sample1     > tmp; mv tmp $se_jcec_ijcend_sample1
paste -d ' ' $se_jcec_idend_sample1 $se_jcec_sjcend_sample1     > tmp; mv tmp $se_jcec_sjcend_sample1
paste -d ' ' $se_jcec_idend_sample1 $se_jcec_inclenend_sample1  > tmp; mv tmp $se_jcec_inclenend_sample1
paste -d ' ' $se_jcec_idend_sample1 $se_jcec_skiplenend_sample1 > tmp; mv tmp $se_jcec_skiplenend_sample1

paste -d ' ' $se_jc_idend_sample2 $se_jc_incend_sample2         > tmp; mv tmp $se_jc_incend_sample2
paste -d ' ' $se_jc_idend_sample2 $se_jc_ijcend_sample2         > tmp; mv tmp $se_jc_ijcend_sample2
paste -d ' ' $se_jc_idend_sample2 $se_jc_sjcend_sample2         > tmp; mv tmp $se_jc_sjcend_sample2
paste -d ' ' $se_jc_idend_sample2 $se_jc_inclenend_sample2      > tmp; mv tmp $se_jc_inclenend_sample2
paste -d ' ' $se_jc_idend_sample2 $se_jc_skiplenend_sample2     > tmp; mv tmp $se_jc_skiplenend_sample2

paste -d ' ' $se_jcec_idend_sample2 $se_jcec_incend_sample2     > tmp; mv tmp $se_jcec_incend_sample2
paste -d ' ' $se_jcec_idend_sample2 $se_jcec_ijcend_sample2     > tmp; mv tmp $se_jcec_ijcend_sample2
paste -d ' ' $se_jcec_idend_sample2 $se_jcec_sjcend_sample2     > tmp; mv tmp $se_jcec_sjcend_sample2
paste -d ' ' $se_jcec_idend_sample2 $se_jcec_inclenend_sample2  > tmp; mv tmp $se_jcec_inclenend_sample2
paste -d ' ' $se_jcec_idend_sample2 $se_jcec_skiplenend_sample2 > tmp; mv tmp $se_jcec_skiplenend_sample2
#------------------------------------------------------------------------------------------------------
#
# A3SS 
#
#------------------------------------------------------------------------------------------------------
#
# these can be used for sashimi plots later
#
#sample1="SAMPLE_1/REP_1/aligned.sorted.bam"
#sample2="SAMPLE_2/REP_1/aligned.sorted.bam"

#cp $dir$sample1 $sample1_aligned_sorted_bam
#cp $dir$sample2 $sample2_aligned_sorted_bam
#------------------------------------------------------------------------------------------------------
#
# fini
#
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
