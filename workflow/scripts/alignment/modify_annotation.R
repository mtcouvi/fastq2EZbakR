#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Modify user-provided annotation to include a full-length gene
## transcript which when quantified will provide information about
## intronic coverage levels to simulate.

# Load dependencies ------------------------------------------------------------

library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(optparse)

# Process parameters -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-g", "--gtf", type = "character"),
              help = "Path to user-provided annotation"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to modified annotation output"),
  make_option(c("-c", "--clean"),
              action = "store_true",
              default = FALSE,
              help = "Remove KI and GL chromosomes")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Modify annotation ------------------------------------------------------------


# Load initial gtf
gtf <- as_tibble(rtracklayer::import(opt$gtf))

# Identify genes with introns
intronless <- gtf %>%
  filter(type == "exon") %>%
  group_by(transcript_id, gene_id) %>%
  dplyr::count() %>%
  filter(n == 1) %>%
  ungroup() %>%
  dplyr::select(gene_id) %>%
  dplyr::distinct()

# Full genes to add to annotation
genes <- gtf %>%
  filter(type == "transcript") %>%
  filter(!(gene_id %in% intronless$gene_id)) %>%
  group_by(gene_id) %>%
  summarise(seqnames = unique(seqnames),
            start = min(start),
            end = max(end),
            strand = unique(strand)) %>%
  mutate(transcript_id = paste0(gene_id, ".I"),
         width = (end - start) + 1,
         source = "NRsim",
         score = NA,
         phase = NA,
         type = "transcript")

# Duplicate to call it an exon
exons <- genes %>%
  mutate(type = "exon")

# Generate final gtf
final_gtf <- bind_rows(list(gtf, genes, exons))


# Remove garbage chromosomes
if(opt$clean){

  final_gtf <- final_gtf %>%
    filter(!grepl("KI", seqnames) & !grepl("GL", seqnames))
    
}


# Save final gtf
final_gr <- GRanges(seqnames = Rle(final_gtf$seqnames),
                    ranges = IRanges(final_gtf$start, end = final_gtf$end),
                    strand = Rle(final_gtf$strand))

mcols(final_gr) <- final_gtf %>%
  dplyr::select(-seqnames, -start, -end, -strand, -width)

rtracklayer::export(final_gr, con = opt$output)
