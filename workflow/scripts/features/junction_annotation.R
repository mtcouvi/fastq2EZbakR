#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Craft a junction feature that featureCounts could
## assign reads to in fastq2EZbakR 

# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(optparse)
library(tidyr)


args = commandArgs(trailingOnly = TRUE)


option_list <- list(
  make_option(c("-r", "--reference", type="character"),
              help = "Path to input gtf file to be expanded"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to output gtf file"),
  make_option(c("-e", "--eeoverhang", type = "double"),
              default = 5,
              help = "Overhang for exon-exon junctions"),
  make_option(c("-i", "--eioverhang", type = "double"),
              default = 5,
              help = "Overhang for exon-intron junctions"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.

overhang <- opt$eeoverhang
ei_overhang <- opt$eioverhang


# Mess around with an annotation -----------------------------------------------

gtf <- rtracklayer::import(opt$reference)


gtfdf <- as_tibble(gtf[mcols(gtf)$type == "exon"]) %>%
  arrange(gene_id, transcript_id, start) %>%
  group_by(gene_id, transcript_id) %>%
  mutate(exon_number = 1:n())


### Create annotation of exon-exon junctions
junctions <- gtfdf %>%
  ungroup() %>%
  dplyr::select(seqnames, start, end, strand, gene_id, transcript_id, exon_number) %>%
  dplyr::distinct() %>%
  group_by(gene_id, transcript_id) %>%
  dplyr::arrange(gene_id, transcript_id, exon_number) %>%
  dplyr::mutate(fiveprime_start = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = end - overhang
  ),
  fiveprime_end = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = end
  ),
  threeprime_end = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = lead(start, order = exon_number) + overhang
  ),
  threeprime_start = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = lead(start, order = exon_number)
  ),
  junction_id = exon_number) %>%
  na.omit() %>%
  dplyr::ungroup() %>%
  dplyr::select(-start, -end, -exon_number, -transcript_id, -junction_id) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(fiveprime_start) %>%
  dplyr::mutate(junction_id = paste0(gene_id, ".J", 1:n())) %>%
  dplyr::ungroup() %>%
  pivot_longer(
    cols = contains("start") | contains("end"),
    names_to = c("type", ".value"),
    names_pattern = "(fiveprime_|threeprime_)(.*)"
  ) %>%
  dplyr::mutate(type = "eej",
                width = overhang - 1)


### Create annotation of exon-intron junctions
eij <- gtfdf %>%
  dplyr::mutate(exon_number = as.numeric(exon_number)) %>%
  dplyr::select(seqnames, start, end, strand, gene_id, transcript_id, exon_number) %>%
  dplyr::distinct() %>%
  group_by(gene_id, transcript_id) %>%
  dplyr::arrange(gene_id, transcript_id, exon_number) %>%
  dplyr::mutate(fiveprime_start = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = end + 1
  ),
  fiveprime_end = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = end + ei_overhang + 1
  ),
  threeprime_end = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = lead(start, order = exon_number) - 1
  ),
  threeprime_start = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = lead(start, order = exon_number) - (ei_overhang + 1)
  )) %>%
  na.omit() %>%
  dplyr::ungroup() %>%
  dplyr::select(-start, -end, -exon_number, -transcript_id) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(fiveprime_start) %>%
  dplyr::mutate(junction_id = paste0(gene_id, ".J", 1:n())) %>%
  dplyr::ungroup() %>%
  pivot_longer(
    cols = contains("start") | contains("end"),
    names_to = c("type", ".value"),
    names_pattern = "(fiveprime_|threeprime_)(.*)"
  ) %>%
  dplyr::mutate(type = "eij",
                width = ei_overhang - 1)


final_gtf <- bind_rows(list(as_tibble(gtf), 
                            junctions,
                            eij))

# Convert to GenomicRanges object and export
new_GR <- GRanges(seqnames = Rle(final_gtf$seqnames),
                  ranges = IRanges(final_gtf$start, end = final_gtf$end, 
                                   names = 1:nrow(final_gtf)),
                  strand = Rle(final_gtf$strand))

mcols(new_GR) <- final_gtf %>%
  dplyr::select(-seqnames, -start, -end, -width, -strand)

rtracklayer::export(new_GR,
                    con = opt$output)
