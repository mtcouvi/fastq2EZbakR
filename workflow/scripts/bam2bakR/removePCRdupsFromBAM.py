#!/usr/bin/env python
"""

Date : March 18, 2016

Author : Heather Landry

Remove reads in bam file resulting from PCR duplicates. This script records all barcodes and coordinates
at a specific position. For every bam line, if the barcode and coordinate has not been seen previously, it
will print; if the barcode and position has been seen previously, it will not print to a new file.
                                        
Update April 25, 2017  - Mary Couvillion:
Skip lines where target id (tid) is <0
Update Sept 2018 - Mary Couvillion:
add CIGAR string to UMI
Update May 2025 - Mary Couvillion:
add flag to UMI to make sure both reads of paired end are kept even if they have the same coordinates and CIGAR                 
"""
import sys, pysam

input_bam = snakemake.input.input_bam
output_bam = snakemake.output.output_bam

iBAM = pysam.Samfile(input_bam, 'rb')
oBAM = pysam.Samfile(output_bam, 'wb', template=iBAM)

MB = set()

# read through starting bam file
for read in iBAM:
    mb = read.qname.split('_MolecularBarcode:')[1]
    if read.tid < 0:
        continue
    if read.tid >= 0:
        chrom = iBAM.getrname(read.tid)
        cigar = read.cigarstring
        flag = read.flag
        
    # selecting the 5' position for pos strand
    if not read.is_reverse:
        start = read.reference_start
        std='pos'
    
    # selecting the 5' position for neg strand
    if  read.is_reverse:
        start = read.reference_end
        std='neg'

    key = str(chrom)+"_"+str(start)+"_"+str(std)+"_"+str(mb)+"_"+str(cigar)+"_"+str(flag)
    
    # output 1 read per molecular barcode
    if key not in MB:
        MB.add(key)
        oBAM.write(read)

iBAM.close()
oBAM.close()

# module load dev/python/2.7.10

