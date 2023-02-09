#!/usr/bin/env Rscript

# print out working directory to make sure you are in the correct place
getwd()

# print out date and time for reference 
Sys.time()
Sys.Date()

# set options and load libraries
options(stringsAsFactors = FALSE)

# load in config and extract data
config <- read.table("config_for_rmats_and_postprocessing.txt")
colnames(config) <- c("parameter", "value")

# get genome version
ref_genome <- basename(config[config$parameter == "fasta", "value"])
genome_version <- strsplit(ref_genome, "[.]")[[1]][1]

if(genome_version == "GRCh38"){
  print("Loading human genome, hg38")
  if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

  library(BSgenome.Hsapiens.UCSC.hg38)
  genome <- BSgenome.Hsapiens.UCSC.hg38
}

if(genome_version == "GRCm38"){
  print("Loading mouse genome, mm10")
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- BSgenome.Mmusculus.UCSC.mm10
}
if(!(genome_version %in% c("GRCh38",  "GRCm38"))){
  print("ERROR! Need to install new genome library. Check BSgenome library")
}
if (!requireNamespace("BiocManager", quietly = TRUE))
#chooseCRANmirror()
#source("https://bioconductor.org/biocLite.R")
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("GenomicRanges", ask = FALSE)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", ask = FALSE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

genome <- BSgenome.Hsapiens.UCSC.hg38

# load in JCEC files for each event type 
CA <- read.table("SE.MATS.JCEC.txt", header=TRUE, sep="\t")
print("CA loaded")
RI <- read.table("RI.MATS.JCEC.txt", header=TRUE, sep="\t")
print("RI loaded")
A3SS <- read.table("A3SS.MATS.JCEC.txt", header=TRUE, sep="\t")
print("A3SS loaded")
A5SS <- read.table("A5SS.MATS.JCEC.txt", header=TRUE, sep="\t")
print("A5SS loaded")
MXE <- read.table("MXE.MATS.JCEC.txt", header=TRUE, sep="\t")
print("MXE loaded")

# get names of b1 and b2
b1 <- basename(config[config$parameter=="b1", "value"])
b1 <- strsplit(b1, "[.]")[[1]][1]
b2 <- basename(config[config$parameter=="b2", "value"])
if(length(b2 > 0)){
b2 <- strsplit(b2, "[.]")[[1]][1]
} else {
b2 <- "none"
print("b2 not present")
}

# get names of gtf used
rmats_gtf <- basename(config[config$parameter == "rmats_gtf", "value"])
orig_gtf <-  basename(config[config$parameter == "ref_gtf", "value"])

# get rmats id
rmats_id <-  basename(config[config$parameter == "rmats_id", "value"])

# function used to calculate mean for Inc and PSI
calc_mean <- function(val){
  new_val <- as.numeric(strsplit(val, ";")[[1]])
  final_val <- mean(new_val, na.rm=TRUE)
  return(final_val)
}
chr_list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",  "chr11",
              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
              "chrX", "chrY", "chrM")

# main script to generate output matrix
get_postprocessing_table <- function(df){
  # get event type column
  event_type <- deparse(substitute(df))
  print(sprintf("Processing %s....", event_type))
  print("Table of number of events per chromosome")
  print(table(df$chr))
  print("NOTE: We removes caffolds, assembly patches and alternate loci (haplotypes). These were needed for mapping purposes.")
  df <- subset(df, chr %in% chr_list)
  df$Event_type <-event_type
  #add 1 to coordinates that are zero base
  if(event_type=="CA"){
    df$exonStart_0base <- (df$exonStart_0base +1)
    df$upstreamES <- (df$upstreamES +1)
    df$downstreamES <- (df$downstreamES +1)
  }
  if(event_type=="RI"){
    df$riExonStart_0base <- (df$riExonStart_0base +1)
    df$upstreamES <- (df$upstreamES +1)
    df$downstreamES <- (df$downstreamES +1)
  }
  if(event_type=="A3SS" | event_type=="A5SS"){
    df$longExonStart_0base <- (df$longExonStart_0base +1)
    df$shortES <- (df$shortES +1)
    df$flankingES <- (df$flankingES +1)
  }
  if(event_type=="MXE"){
    df$X1stExonStart_0base <- (df$X1stExonStart_0base +1)
    df$X2ndExonStart_0base <- (df$X2ndExonStart_0base +1)
    df$upstreamES <- (df$upstreamES +1)
    df$downstreamES <- (df$downstreamES +1)
  }
  # add long and short IDs, and get exon coordinates
  if(event_type=="CA"){
    # create long and short IDs
    df$long_ID <- ifelse(df$strand == "+", 
                         paste(df$Event_type, df$geneSymbol, df$chr, df$strand,df$upstreamES, df$upstreamEE,df$exonStart_0base, df$exonEnd, df$downstreamES, df$downstreamEE, sep=";" ),
                         paste(df$Event_type, df$geneSymbol, df$chr, df$strand, df$downstreamEE, df$downstreamES,df$exonEnd,df$exonStart_0base, df$upstreamEE, df$upstreamES, sep=";"))
    df$short_ID <- ifelse(df$strand == "+",
                          paste(df$Event_type, df$geneSymbol, df$chr, df$strand, df$upstreamEE,df$exonStart_0base, df$exonEnd, df$downstreamES, sep=";"),
                          paste(df$Event_type, df$geneSymbol, df$chr, df$strand, df$downstreamES,df$exonEnd,df$exonStart_0base, df$upstreamEE, sep=";"))
    # get coordinates, differnt for - strand
    df$E1_start <- ifelse(df$strand == "+", df$upstreamES, df$downstreamES)
    df$E1_end <- ifelse(df$strand == "+", df$upstreamEE, df$downstreamEE)
    df$E2_start <- df$exonStart_0base
    df$E2_end <- df$exonEnd
    df$E3_start <- ifelse(df$strand == "+", df$downstreamES, df$upstreamES)
    df$E3_end <- ifelse(df$strand == "+", df$downstreamEE, df$upstreamEE)
  }
  if(event_type=="RI"){
    # make long and short IDs
    df$long_ID <- ifelse(df$strand == "+",
                    paste(df$Event_type, df$geneSymbol, df$chr, df$strand, df$upstreamES, df$upstreamEE, df$upstreamEE+1, df$downstreamES-1, df$downstreamES, df$downstreamEE, sep=";" ),
                    paste(df$Event_type, df$geneSymbol, df$chr, df$strand, df$downstreamEE, df$downstreamES, df$downstreamES-1, df$upstreamEE+1, df$upstreamEE, df$upstreamES, sep=";" ))
    df$short_ID <- ifelse(df$strand == "+", 
                          paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$upstreamEE, df$upstreamEE+1, df$downstreamES-1, df$downstreamES, sep=";"),
                          paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$downstreamES, df$downstreamES-1, df$upstreamEE+1, df$upstreamEE, sep=";"))
    # get coordinates, different for - strand 
    df$E1_start <- ifelse(df$strand == "+", df$upstreamES, df$downstreamES)
    df$E1_end <- ifelse(df$strand == "+", df$upstreamEE, df$downstreamEE)
    df$E2_start <- df$upstreamEE+1
    df$E2_end <- df$downstreamES-1
    df$E3_start <- ifelse(df$strand == "+", df$downstreamES, df$upstreamES)
    df$E3_end <- ifelse(df$strand == "+", df$downstreamEE, df$upstreamEE)
  }
  if(event_type=="A3SS"){ 
    # make long and short IDs
    df$long_ID <- ifelse(df$strand == "+", 
                  paste(df$Event_type, df$geneSymbol,df$chr, df$strand, df$flankingES, df$flankingEE, df$longExonStart_0base, df$shortES-1,df$shortES, df$shortEE, sep=";" ),
                  paste(df$Event_type, df$geneSymbol,df$chr, df$strand, df$flankingEE, df$flankingES, df$longExonEnd, df$shortEE+1, df$shortEE, df$shortES, sep=";" ))
    
    df$short_ID <- ifelse(df$strand == "+", 
                          paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$flankingEE, df$longExonStart_0base, df$shortES-1,df$shortES, sep=";"),
                          paste(df$Event_type,df$geneSymbol, df$chr, df$strand,  df$flankingES, df$longExonEnd, df$shortEE+1, df$shortEE, sep=";"))
    # get coordinates, different for - strand 
    df$E1_start <- df$flankingES
    df$E1_end <- df$flankingEE
    df$E2_start <- ifelse(df$strand == "+", df$longExonStart_0base, df$shortEE+1)
    df$E2_end <- ifelse(df$strand == "+", df$shortES-1, df$longExonEnd)
    df$E3_start <-  df$shortES
    df$E3_end <- df$shortEE
  }
  if(event_type=="A5SS"){
    # make long and short IDs
    df$long_ID <- ifelse(df$strand == "+", 
                         paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$shortES, df$shortEE, df$shortEE+1, df$longExonEnd, df$flankingES, df$flankingEE, sep=";"),
                         paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$shortEE, df$shortES, df$shortES-1, df$longExonStart_0base, df$flankingEE, df$flankingES, sep=";"))
    df$short_ID <- ifelse(df$strand == "+", 
                          paste(df$Event_type, df$geneSymbol,df$chr, df$strand, df$shortEE, df$shortEE+1, df$longExonEnd, df$flankingES, sep=";"),
                          paste(df$Event_type, df$geneSymbol,df$chr, df$strand, df$shortES, df$shortES-1, df$longExonStart_0base, df$flankingEE, sep=";"))
    # get coordinates, diff for - strand
    df$E1_start <- df$shortES
    df$E1_end <- df$shortEE
    df$E2_start <- ifelse(df$strand == "+", df$shortEE+1, df$longExonStart_0base)
    df$E2_end <- ifelse(df$strand == "+", df$longExonEnd, df$shortES-1)
    df$E3_start <- df$flankingES
    df$E3_end <- df$flankingEE
  }
  if(event_type =="MXE"){
    df$long_ID <- ifelse(df$strand == "+",
                         paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$upstreamES, df$upstreamEE, df$X1stExonStart_0base, df$X1stExonEnd,
                               df$X2ndExonStart_0base, df$X2ndExonEnd,  df$downstreamES, df$downstreamEE, sep=";"), 
                         paste(df$Event_type,df$geneSymbol, df$chr, df$strand, df$downstreamEE, df$downstreamES, df$X2ndExonEnd,df$X2ndExonStart_0base,  
                              df$X1stExonEnd,df$X1stExonStart_0base, df$upstreamEE, df$upstreamES, sep=";"))
    df$short_ID <- ifelse(df$strand == "+", 
                          paste(df$Event_type, df$geneSymbol,df$chr, df$strand,df$upstreamEE, df$X1stExonStart_0base, df$X1stExonEnd,
                                df$X2ndExonStart_0base, df$X2ndExonEnd,  df$downstreamES ,  sep=";"), 
                          paste(df$Event_type, df$geneSymbol,df$chr, df$strand, df$downstreamES, df$X2ndExonEnd,df$X2ndExonStart_0base,  
                                df$X1stExonEnd,df$X1stExonStart_0base, df$upstreamEE,  sep=";"))
    # get coordinates, differnt for - strand
    df$E1_start <- ifelse(df$strand == "+", df$upstreamES, df$downstreamES)
    df$E1_end <- ifelse(df$strand == "+", df$upstreamEE, df$downstreamEE)
    df$E2_start <- ifelse(df$strand == "+", df$X1stExonStart_0base, df$X2ndExonStart_0base)
    df$E2_end <- ifelse(df$strand == "+", df$X1stExonEnd, df$X2ndExonEnd)
    df$E3_start <- ifelse(df$strand == "+", df$X2ndExonStart_0base, df$X1stExonStart_0base)
    df$E3_end <- ifelse(df$strand == "+", df$X2ndExonEnd, df$X1stExonEnd)
    df$E4_start <- ifelse(df$strand == "+", df$downstreamES, df$upstreamES)
    df$E4_end <- ifelse(df$strand == "+", df$downstreamEE, df$upstreamEE)
  }
  # get sequences 
  if(event_type != "MXE"){
    # get coordinates
    df$E1 <- paste0(df$chr, ":", df$E1_start, "-", df$E1_end)
    df$E2 <- paste0(df$chr, ":", df$E2_start, "-", df$E2_end)
    df$E3 <- paste0(df$chr, ":", df$E3_start, "-", df$E3_end)
    df$E4 <- NA
    # make ranges object
    E1_ranges <- GRanges(df$chr, IRanges(start = df$E1_start, end = df$E1_end), strand = df$strand)
    E2_ranges <- GRanges(df$chr, IRanges(start = df$E2_start, end = df$E2_end), strand = df$strand)
    E3_ranges <- GRanges(df$chr, IRanges(start = df$E3_start, end = df$E3_end), strand = df$strand)
    # get sequences
    E1_seq <- getSeq(genome, E1_ranges)
    df$S1 <- as.character(E1_seq)
    E2_seq <- getSeq(genome, E2_ranges)
    df$S2 <- as.character(E2_seq)
    E3_seq <- getSeq(genome, E3_ranges)
    df$S3 <- as.character(E3_seq)
    df$S4 <- NA
  }
  else{
    # get coordinates
    df$E1 <- paste0(df$chr, ":", df$E1_start, "-", df$E1_end)
    df$E2 <- paste0(df$chr, ":", df$E2_start, "-", df$E2_end)
    df$E3 <- paste0(df$chr, ":", df$E3_start, "-", df$E3_end)
    df$E4 <- paste0(df$chr, ":", df$E4_start, "-", df$E4_end)
    # make ranges object
    E1_ranges <- GRanges(df$chr, IRanges(start = df$E1_start, end = df$E1_end), strand = df$strand)
    E2_ranges <- GRanges(df$chr, IRanges(start = df$E2_start, end = df$E2_end), strand = df$strand)
    E3_ranges <- GRanges(df$chr, IRanges(start = df$E3_start, end = df$E3_end), strand = df$strand)
    E4_ranges <- GRanges(df$chr, IRanges(start = df$E4_start, end = df$E4_end), strand = df$strand)
    # get sequences
    E1_seq <- getSeq(genome, E1_ranges)
    df$S1 <- as.character(E1_seq)
    E2_seq <- getSeq(genome, E2_ranges)
    df$S2 <- as.character(E2_seq)
    E3_seq <- getSeq(genome, E3_ranges)
    df$S3 <- as.character(E3_seq)
    E4_seq <- getSeq(genome, E4_ranges)
    df$S4 <- as.character(E4_seq)
  }
  # create blat column
  df$blat <- paste(">", df$short_ID, sep = "")
  df$blat <- paste(df$blat, df$S1, sep="\r")
  df$blat <- paste(df$blat, df$S2, sep="\r")
  df$blat <- paste(df$blat, df$S3, sep="\r")
  if(event_type =="MXE"){
    df$blat <- paste(df$blat, df$S4, sep="\r")
  }
  # replace ',' with ';'
  df[,c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel1", "IncLevel2")] <- apply(df[,c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel1", "IncLevel2")], MARGIN=c(1,2), function(x) {
    gsub(",", ";", x)
  })
  # calculate mean values
  df$mIC1 <-  sapply(df$IJC_SAMPLE_1, function(x) calc_mean(x))
  df$mSC1 <- sapply(df$SJC_SAMPLE_1, function(x) calc_mean(x))
  df$mIC2 <- sapply(df$IJC_SAMPLE_2, function(x) calc_mean(x))
  df$mSC2 <- sapply(df$SJC_SAMPLE_2, function(x) calc_mean(x))
  df$mPSI1 <- sapply(df$IncLevel1, function(x) calc_mean(x))
  df$mPSI2 <- sapply(df$IncLevel2, function(x) calc_mean(x))
  df$IncLevelDifference <- df$IncLevelDifference*(-1)
  # reorder and rename columns, keeping only the ones we want
  final_df <- df[,c("Event_type", "short_ID", "long_ID", "geneSymbol", "GeneID", "strand", "IJC_SAMPLE_1", "SJC_SAMPLE_1", 
                    "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel1", "IncLevel2", "mIC1", "mSC1", "mIC2", "mSC2", "mPSI1", "mPSI2",
                    "IncLevelDifference", "PValue", "FDR", "E1", "E2", "E3", "E4", "S1", "S2", "S3", "S4", "blat", "IncFormLen", "SkipFormLen")]
  colnames(final_df) <- c("Event_type", "short_ID", "long_ID", "geneSymbol", "geneID", "strand", "IC1", "SC1", 
                          "IC2", "SC2", "PSI1", "PSI2", "mIC1", "mSC1", "mIC2", "mSC2", "mPSI1", "mPSI2",
                          "dPSI_2_minus_1", "PValue", "FDR", "E1", "E2", "E3", "E4", "S1", "S2", "S3", "S4", "blat", "IncFormLen", "SkipFormLen")
  print(sprintf("%s done.", event_type))
  return(final_df)
}
CA <- get_postprocessing_table(CA)
RI <- get_postprocessing_table(RI)
A3SS <- get_postprocessing_table(A3SS)
A5SS <- get_postprocessing_table(A5SS)
MXE <- get_postprocessing_table(MXE)

# merge all together
all_events <- rbind(CA, RI, A3SS, A5SS, MXE)
all_events$genome <- ref_genome
all_events$rmats_gtf <- rmats_gtf
all_events$orig_gtf <- orig_gtf

# write out file
date <- Sys.Date()
sample_prefix <- rmats_id
if(rmats_gtf=="gffcmp.annotated.corrected.gtf"){
  write.csv(all_events, paste0("Unfiltered_","mergedGTF_", sample_prefix, "_", date, ".csv"),row.names=FALSE)
}
if(rmats_gtf!="gffcmp.annotated.corrected.gtf"){
  write.csv(all_events, paste0("Unfiltered_", sample_prefix, "_", date, ".csv"),row.names=FALSE)
}

# filter events
filter_events <- function(df){
  print("Filtering")
  print(sprintf("There were %i events in the unfiltered file", nrow(df)))
  new_df <- subset(df, FDR <= 0.05)
  new_df <- subset(new_df, abs(dPSI_2_minus_1) >= 0.1)
  new_df <- subset(new_df, (mIC1 >= 5 | mIC2 >=5) & (mSC1 >=5 | mSC2 >= 5))
  print(sprintf("There were %i events in the filtered file", nrow(new_df)))
  return(new_df)
}
all_filt <- filter_events(all_events)
if(rmats_gtf=="gffcmp.annotated.corrected.gtf"){
  write.csv(all_filt, paste0("Filtered_","mergedGTF_", sample_prefix, "_", date, ".csv"),row.names=FALSE)
}
if(rmats_gtf!="gffcmp.annotated.corrected.gtf"){
  write.csv(all_filt, paste0("Filtered_", sample_prefix, "_", date, ".csv"),row.names=FALSE)
}

if(nrow(all_filt)>100){                     
# get summary table
print("Generating summary table")
get_summary_table <- function(unfilt_df, filt_df){
  # get unfiltered events
  unfilt_table <- as.data.frame(table(unfilt_df$Event_type))
  colnames(unfilt_table) <- c("Event_type", "Detected")
  # filtered events
  filt_table <- as.data.frame(table(filt_df$Event_type))
  colnames(filt_table) <- c("Event_type", "Filtered")
  # pos and neg PSI tables
  pos_table <- as.data.frame(table(subset(filt_df, dPSI_2_minus_1 > 0)$Event_type))
  colnames(pos_table) <- c("Event_type", "pos_PSI")
  neg_table <- as.data.frame(table(subset(filt_df, dPSI_2_minus_1 < 0)$Event_type))
  colnames(neg_table) <- c("Event_type", "neg_PSI")
  # merge together
  final_df <- Reduce(function(x,y) merge(x=x, y=y, by="Event_type", all=TRUE), 
                     list(unfilt_table, filt_table, pos_table, neg_table))
  return(final_df)
}
summary_table <- get_summary_table(all_events, all_filt)
if(rmats_gtf=="gffcmp.annotated.corrected.gtf"){
  write.csv(summary_table, paste0("SummaryTable_","mergedGTF_", sample_prefix, "_", date, ".csv", sep=""),row.names=FALSE)
}
if(rmats_gtf!="gffcmp.annotated.corrected.gtf"){
  write.csv(summary_table, paste0("SummaryTable_", sample_prefix, "_", date, ".csv", sep=""),row.names=FALSE)
}
}
                     
# generate bam file list
b1_file <- config[config$parameter=="b1", "value"]
b2_file <- config[config$parameter=="b2", "value"]

b1_list <- t(read.table(b1_file, sep=","))
if(b2 != "none"){
b2_list <- t(read.table(b2_file, sep=","))
bam_list_df <- merge(b1_list, b2_list, by="row.names", suffixes = c("_b1", "_b2"), all=TRUE)
bam_list_df <- bam_list_df[,-1]
colnames(bam_list_df) <- lapply(colnames(bam_list_df), function(x){
  new_val <- strsplit(x, "_")[[1]][2]
  return(new_val)
})
} else {
bam_list_df <- b1_list
}

if(rmats_gtf=="gffcmp.annotated.corrected.gtf"){
  write.csv(bam_list_df, paste0("BamList_","mergedGTF_", sample_prefix, "_", date, ".csv", sep=""),row.names=FALSE)
}
if(rmats_gtf!="gffcmp.annotated.corrected.gtf"){
  write.csv(bam_list_df, paste0("BamList_", sample_prefix, "_", date, ".csv", sep=""),row.names=FALSE)
}

