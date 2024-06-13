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
  make_option(c("--starjunc", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to junctions via relevant STAR tags."),              
  make_option(c("-j", "--eej", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exon-exon junctions"),
  make_option(c("--eij", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exon-intron junctions"),
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

sample <- paste0("^", opt$sample, "_counts")

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
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
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
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
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
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
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
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
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

  sample <- paste0("^", opt$sample, ".csv")
  
  transcripts_file <- list.files("./results/read_to_transcripts/",
                              pattern = sample, full.names = TRUE)[1]
  
  transcripts <- fread(transcripts_file)
  
  colnames(transcripts) <- c("qname", "bamfile_transcripts")
  
  
  transcripts <- transcripts[ , c("qname", "bamfile_transcripts")]  
  
  setkey(transcripts, qname)

  muts <- transcripts[muts]
  
  feature_vect <- c(feature_vect, "bamfile_transcripts")

}


# Exon-exon junction assignment
if(opt$eej){
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
  transcripts_file <- list.files("./results/featurecounts_eej/",
                              pattern = sample, full.names = TRUE)[1]
  
  transcripts <- fread(transcripts_file)
  
  colnames(transcripts) <- c("qname", "status", "nhits", "ee_junction_id")
  
  
  transcripts <- transcripts[ nhits > 0 , c("qname", "ee_junction_id")]
  
  transcripts[, ee_junction_id := gsub(",", "+", ee_junction_id)]
  
  setkey(transcripts, qname)
  
  muts <- transcripts[muts]
  
  feature_vect <- c(feature_vect, "ee_junction_id")

  
}


# Exon-intron junction assignment
if(opt$eij){
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
  transcripts_file <- list.files("./results/featurecounts_eij/",
                              pattern = sample, full.names = TRUE)[1]
  
  transcripts <- fread(transcripts_file)
  
  colnames(transcripts) <- c("qname", "status", "nhits", "ei_junction_id")
  
  
  transcripts <- transcripts[ nhits > 0 , c("qname", "ei_junction_id")]
  
  transcripts[, ei_junction_id := gsub(",", "+", ei_junction_id)]
  
  setkey(transcripts, qname)
  
  muts <- transcripts[muts]
  
  feature_vect <- c(feature_vect, "ei_junction_id")

  
}


# STAR junction assignment
if(opt$starjunc){
  
  sample <- paste0("^", opt$sample, ".csv.gz")
  
  transcripts_file <- list.files("./results/read_to_junctions/",
                              pattern = sample, full.names = TRUE)[1]
  

  message("File is:")
  print(transcripts_file)

  message("file looks like:")

  transcripts <- fread(transcripts_file)
  
  print(head(transcripts))

  colnames(transcripts) <- c("qname", "junction_start", "junction_end")
    
  message("file looks like:")


  setkey(transcripts, qname)
  
  print(head(transcripts))


  muts <- transcripts[muts]
  
  feature_vect <- c(feature_vect, "junction_start", "junction_end")

}


# Write to final output
write_csv(muts,
          file = opt$output)


##### MAKE CB

muts_to_keep <- unlist(strsplit(opt$muttypes, ","))
bases_to_keep <- paste0("n", substr(muts_to_keep, start = 1, stop = 1))

cols_to_keep <- c("sample", "rname", feature_vect, muts_to_keep, bases_to_keep, "sj")

muts[, sample := opt$sample]

print("cols_to_keep looks like:")
print(cols_to_keep)

muts <- muts[, .(n = .N), by = cols_to_keep]


write_csv(muts,
          file = opt$cBoutput)
