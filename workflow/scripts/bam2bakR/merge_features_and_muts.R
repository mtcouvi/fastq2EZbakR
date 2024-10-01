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
              help = "Sample name"),
  make_option(c("--annotation", type = "character"),
              help = "Path to annotation GTF file.")          
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Combine tables --------------------------------------------------------------

sample <- paste0("^", opt$sample, "_counts")

print(paste0("sample is: ", opt$sample))


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

  muts[, GF := ifelse(is.na(GF), "__no_feature", GF)]
  
  feature_vect <- c(feature_vect, "GF")
  
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

  muts[, bamfile_transcripts := ifelse(is.na(bamfile_transcripts), "__no_feature", bamfile_transcripts)]


  if(opt$genes){
    ### Filter out incorrectly assigned isoforms
    ### due to STAR assigning read to isoform on oppposite
    ### strand.
    ### Idea is that featureCount's gene assignment correctly
    ### assigns the gene, so the annotation can be used to determine
    ### the set of isoforms that are legit.
    library(rtracklayer)
    library(tidyr)
    library(dplyr)

    # Table of set of isoforms from each gene
    gene2transcript <- rtracklayer::import(opt$annotation) %>% 
      dplyr::as_tibble() %>%
      dplyr::filter(type == "transcript") %>%
      dplyr::select(gene_id, transcript_id) %>%
      dplyr::distinct() %>%
      dplyr::rename(GF = gene_id)

    setDT(gene2transcript)
    setkey(gene2transcript, GF, transcript_id)

    # Unique set of bamfile_transcripts and GF combos
    current_assignments <- muts %>%
      dplyr::select(GF, bamfile_transcripts) %>%
      dplyr::distinct() %>%
      dplyr::mutate(transcript_id = bamfile_transcripts) %>% 
      tidyr::separate_rows(transcript_id, sep = "\\+")


    setDT(current_assignments)
    setkey(current_assignments, GF, transcript_id)

    current_assignments <- current_assignments[gene2transcript, nomatch = NULL] %>%
      dplyr::group_by(GF, bamfile_transcripts) %>%
      dplyr::summarise(new_bft = paste(transcript_id, collapse="+"))

    if(nrow(current_assignments) == 0){
      stop("Something went wrong, current_assignments is empty!")
    }
    
    setDT(current_assignments)
    setkey(current_assignments, GF, bamfile_transcripts)
    setkey(muts, GF, bamfile_transcripts)

    muts <- muts[current_assignments, nomatch = NULL]

    if(nrow(current_assignments) == 0){
      stop("Something went wrong, muts is empty!")
    }

    muts[, bamfile_transcripts := new_bft]
    muts[, new_bft := NULL]

  }

  
  feature_vect <- c(feature_vect, "bamfile_transcripts")

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

  muts[, XF := ifelse(is.na(XF), "__no_feature", XF)]

  
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

  muts[, exon_bin := ifelse(is.na(exon_bin), "__no_feature", exon_bin)]

  
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

  muts[, transcripts := ifelse(is.na(transcripts), "__no_feature", transcripts)]

  
  feature_vect <- c(feature_vect, "transcripts")

  
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
  
  muts[, ee_junction_id := ifelse(is.na(ee_junction_id), "__no_feature", ee_junction_id)]

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
  
  muts[, ei_junction_id := ifelse(is.na(ei_junction_id), "__no_feature", ei_junction_id)]

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
  
  muts[, junction_start := ifelse(is.na(junction_start), "__no_feature", junction_start)]
  muts[, junction_end := ifelse(is.na(junction_end), "__no_feature", junction_end)]


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
