#!/usr/bin/env python
"""


Author : chatGPT and M. Couvillion

Remove reads in bam file resulting from PCR duplicates and keep the one with fewest mismatches (which are presumably sequencing errors
"""
import sys, pysam

input_bam = snakemake.input.input_bam
output_bam = snakemake.output.output_bam

iBAM = pysam.AlignmentFile(input_bam, 'rb')
oBAM = pysam.AlignmentFile(output_bam, 'wb', template=iBAM)

pysam.index(input_bam)

MB = dict()


def flush_MB(MB, oBAM):
    for _, (_, read) in MB.items():
        oBAM.write(read)


for read in iBAM.fetch(until_eof=True):
    if read.is_unmapped:
        continue

    try:
        mb = read.query_name.split('_MolecularBarcode:')[1]
    except IndexError:
        continue

    chrom = iBAM.get_reference_name(read.reference_id)

    # Flush cache when chromosome changes
    if chrom != current_chrom:
        flush_MB(MB, oBAM)
        MB.clear()
        current_chrom = chrom

    # Strand + position
    if not read.is_reverse:
        start = read.reference_start
        std = 'pos'
    else:
        start = read.reference_end
        std = 'neg'

    cigar = read.cigarstring
    flag = read.flag

    key = f"{chrom}_{start}_{std}_{mb}_{cigar}_{flag}"
    try:
        value = read.get_tag("NM")
    except KeyError:
        continue

    if key not in MB or value < MB[key][0]:
        MB[key] = (value, read)

    # Optional: Flush reads if MB gets too big
    if len(MB) > 1_000_000:
        flush_MB(MB, oBAM)
        MB.clear()

# Final flush
flush_MB(MB, oBAM)
iBAM.close()
oBAM.close()

        
# module load dev/python/2.7.10

