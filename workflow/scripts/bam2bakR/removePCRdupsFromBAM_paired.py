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
use info from both reads for paired end data
"""
import sys, pysam

input_bam = snakemake.input.input_bam
output_bam = snakemake.output.output_bam

iBAM = pysam.Samfile(input_bam, 'rb')
oBAM = pysam.Samfile(output_bam, 'wb', template=iBAM)

MB = set()
read_pairs = {}

# read through starting bam file
for read in iBAM:
    # skip unpaired, unmapped, or unmated reads
    if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
        continue

    qname = read.query_name

    # buffer until both reads are available
    if qname not in read_pairs:
        read_pairs[qname] = read
    else:
        mate = read_pairs.pop(qname)

        # Order the reads consistently
        read1, read2 = (read, mate) if not read.is_read2 else (mate, read)

        # Extract molecular barcode
        try:
            mb = read1.query_name.split('_MolecularBarcode:')[1]
        except IndexError:
            continue

        # Get reference name
        chrom1 = iBAM.get_reference_name(read1.reference_id)
        chrom2 = iBAM.get_reference_name(read2.reference_id)

        # Strand-aware 5' start positions
        if not read1.is_reverse:
            start1 = read1.reference_start
        else:
            start1 = read1.reference_end

        if not read2.is_reverse:
            start2 = read2.reference_start
        else:
            start2 = read2.reference_end

        # Create a key using info from both mates
        key = str(chrom1)+"_"+str(start1)+"_"+str(read1.cigarstring)+"_"+str(chrom2)+"_"+str(start2)+"_"+str(read2.cigarstring)+"_"+str(mb)
        
        # Deduplicate based on the full key
        if key not in MB:
            MB.add(key)
            oBAM.write(read1)
            oBAM.write(read2)






    mb = read.qname.split('_MolecularBarcode:')[1]
    if read.tid < 0:
        continue
    if read.tid >= 0:
        chrom = iBAM.getrname(read.tid)
        cigar = read.cigarstring
        
    # selecting the 5' position for pos strand
    if not read.is_reverse:
        start = read.reference_start
        std='pos'
    
    # selecting the 5' position for neg strand
    if  read.is_reverse:
        start = read.reference_end
        std='neg'

    key = str(chrom)+"_"+str(start)+"_"+str(std)+"_"+str(mb)+"_"+cigar
    
    # output 1 read per molecular barcode
    if key not in MB:
        MB.add(key)
        oBAM.write(read)

iBAM.close()
oBAM.close()

# module load dev/python/2.7.10

