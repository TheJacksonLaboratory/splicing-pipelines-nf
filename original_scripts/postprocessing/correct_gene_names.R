#!/usr/bin/env Rscript

library(rtracklayer)

# print out working directory to make sure you are in the correct place
getwd()

# load in new gffcompare gtf
gffcmp_gtf <- import("gffcmp.annotated.gtf", format = "gtf")

# convert to dataframe 
gffcmp_df <- as.data.frame(gffcmp_gtf)

# print out table with class codes
print("This is the summary of code output from gffcompare")
print(table(gffcmp_df$class_code))

# create key dataframe that only has transcripts
key_df <- subset(gffcmp_df, type == "transcript")
print(sprintf("Total transcripts: %i", nrow(key_df)))
print(sprintf("Total with MSTRG names: %i", sum(grepl("MSTRG." ,key_df$gene_id))))
print(sprintf("Total left without gene names after gff compare: %i", sum(is.na(key_df$ref_gene_id))))

# fill in rows where gffcmp did have a gene 
key_df$gene_id[!is.na(key_df$ref_gene_id)] <- key_df$ref_gene_id[!is.na(key_df$ref_gene_id)]

# extract only columns of interest
key_df <- key_df[,c("transcript_id", "gene_id", "gene_name")]

# merge the gffcmp gtf with key to get gene ids for exons 
final_df <- merge(gffcmp_df, key_df, by="transcript_id", all=TRUE, suffixes = c("", "_new"))

# replace old gene name and gene id with new 
final_df$gene_id <- final_df$gene_id_new
final_df$gene_name <- final_df$gene_name_new

#remove extra columns
final_df <- final_df[,1:(ncol(final_df)-2)]

### need to export as gff2
export(final_df, "gffcmp.annotated.corrected.gff", format="gff2")
