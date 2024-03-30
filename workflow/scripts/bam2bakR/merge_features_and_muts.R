#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Merge feature assignment and mutation counts tables

# Load dependencies ------------------------------------------------------------

library(optparse)
library(data.table)
library(readr)

# Process parameters -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-g", "--genes", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to genes"),
  make_option(c("-e", "--exons", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exons"),
  make_option(c("-b", "--exonbins", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exonbins"),
  make_option(c("-t", "--transcripts", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to transcripts"),
  make_option(c("--frombam", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to transcripts from the 
              transcriptome aligned bam file directly. This is more accurate
              than featureCounts based transcript isoform assignment as
              the latter does not account for the splice junctions a read
              is mapped across."),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to full mutation counts/feature assignment output."),
  make_option(c("-c", "--cBoutput", type = "character"),
              help = "Path to cB output; same as full output with some columns averaged out."),
  make_option(c("-m", "--muttypes", type = "character"),
              help = "String of comma separated mutation types to keep in cBs."),
  make_option(c("-s", "--sample", type = "character"),
              help = "Sample name")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Combine tables --------------------------------------------------------------

sample <- paste0(opt$sample, "_counts")

print(paste0("sample is: ", sample))


### Load mutation counts
muts_file <- list.files(path = "./results/counts/",
                        pattern = sample,
                        full.names = TRUE)[1]

muts <- fread(muts_file)

setkey(muts, qname)

feature_vect <- c()

# merge with gene assignments
if(opt$genes){
  
  sample <- paste0(opt$sample, ".s.bam.featureCounts")
  
  genes_file <- list.files("./results/featurecounts_genes/",
                           pattern = sample, full.names = TRUE)[1]
  
  genes <- fread(genes_file)
  
  colnames(genes) <- c("qname", "status", "nhits", "GF")
  
  genes <- genes[ nhits > 0 , c("qname", "GF")]
  
  genes[, GF := gsub(",", "+", GF)]

  setkey(genes, qname)
  
  muts <- genes[muts]
  
  feature_vect <- c(feature_vect, "GF")
  
}


# merge with exon assignments
if(opt$exons){
  
  sample <- paste0(opt$sample, ".s.bam.featureCounts")
  
  exons_file <- list.files("./results/featurecounts_exons/",
                           pattern = sample, full.names = TRUE)[1]
  
  exons <- fread(exons_file)
  
  colnames(exons) <- c("qname", "status", "nhits", "XF")
  
  exons <- exons[ nhits > 0 , c("qname", "XF")]
  
  exons[, XF := gsub(",", "+", XF)]
  
  setkey(exons, qname)

  muts <- exons[muts]
  
  feature_vect <- c(feature_vect, "XF")

}


# Merge with exonbin assignments
if(opt$exonbins){
  
  sample <- paste0(opt$sample, ".s.bam.featureCounts")
  
  exonbins_file <- list.files("./results/featurecounts_exonbins/",
                           pattern = sample, full.names = TRUE)[1]
  
  exonbins <- fread(exonbins_file)
  
  colnames(exonbins) <- c("qname", "status", "nhits", "exon_bin")
  
  
  exonbins <- exonbins[ nhits > 0 , c("qname", "exon_bin")]
  
  
  exonbins[, exon_bin := gsub(",", "+", exon_bin)]
  
  setkey(exonbins, qname)

  muts <- exonbins[muts]
  
  feature_vect <- c(feature_vect, "exon_bin")

}


# Merge with transcript assignments
if(opt$transcripts){
  
  sample <- paste0(opt$sample, ".s.bam.featureCounts")
  
  transcripts_file <- list.files("./results/featurecounts_transcripts/",
                              pattern = sample, full.names = TRUE)[1]
  
  transcripts <- fread(transcripts_file)
  
  colnames(transcripts) <- c("qname", "status", "nhits", "transcripts")
  
  
  transcripts <- transcripts[ nhits > 0 , c("qname", "transcripts")]
  
  transcripts[, transcripts := gsub(",", "+", transcripts)]
  
  setkey(transcripts, qname)
  
  muts <- transcripts[muts]
  
  feature_vect <- c(feature_vect, "transcripts")

  
}

if(opt$frombam){

  sample <- paste0(opt$sample, ".csv")
  
  transcripts_file <- list.files("./results/read_to_transcripts/",
                              pattern = sample, full.names = TRUE)[1]
  
  transcripts <- fread(transcripts_file)
  
  colnames(transcripts) <- c("qname", "bamfile_transcripts")
  
  
  transcripts <- transcripts[ , c("qname", "bamfile_transcripts")]  
  
  setkey(transcripts, qname)

  muts <- transcripts[muts]
  
  feature_vect <- c(feature_vect, "bamfile_transcripts")

}


# Write to final output
write_csv(muts,
          file = opt$output)


##### MAKE CB

muts_to_keep <- strsplit(opt$muttypes, ",")
bases_to_keep <- paste0("n", substr(muts_to_keep, start = 1, stop = 1))

cols_to_keep <- c("sample", feature_vect, muts_to_keep, bases_to_keep)

muts[, sample := opt$sample]

print(paste0("muts_to_keep is: ", muts_to_keep))
print(paste0("bases_to_keep is: ", bases_to_keep))
print(paste0("feature_vect is: ", feature_vect))
print(paste0("cols_to_keep is: ", cols_to_keep))
print("muts looks like:")
head(muts)
key(muts)
setDT(muts)

muts <- muts[, .(n = sum(n)), by = cols_to_keep]


write_csv(muts,
          file = opt$cBoutput)