#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4 # Main bam file that actually matters

# Make sure bam file is query sorted rather than coordinate sorted
    # This is important for mutation counting and other pipeline steps
    # that are parallelized via splitting up the bam file into chunks.
    # Query sorting is required by the custom awk script that splits up 
    # the bam file because the read pairs need to be kept together in the
    # bam file chunks.
# Also filtering out unmapped or unpaired reads
samtools sort -@ "$cpus" -n "$input" | \
samtools fixmate -@ "$cpus" - - | \
samtools view -@ "$cpus" -b -h -q 2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x400 -F 0x800 - > "$output"

echo "* Reads filtered for sample $sample"