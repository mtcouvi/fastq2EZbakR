#!/usr/bin/env python
"""


Author : chatGPT and M. Couvillion

Remove reads in bam file resulting from PCR duplicates and keep, per molecular
barcode family, the read PAIR with the highest combined mean base quality
(presumably the copy least affected by sequencing errors).

Fragment-aware: reads are buffered and paired by query_name, the duplicate
family key is built from BOTH mates (chrom/5'-start/CIGAR of read1 and read2 +
UMI), and both mates of the winning pair are written together so downstream
paired-end mutation calling keeps intact mates. Kept reads are emitted
unmodified, so their MD/NM tags stay valid.

Memory: when the input BAM is coordinate-sorted (e.g. STAR's forced
--outSAMtype BAM SortedByCoordinate), the dictionary of winning pairs is flushed
at every chromosome boundary. This is safe because all PCR duplicates of a
family share identical coordinates and therefore complete (second mate reached)
contiguously, so no further family member can appear once we leave a chromosome.
The half-seen-mate buffer is NEVER flushed at a boundary, since a mate may lie on
a later chromosome. If the input is not coordinate-sorted, flushing is disabled
and all winning pairs are held until the end (always correct, just more memory).
"""
import sys, pysam

input_bam = snakemake.input.input_bam
output_bam = snakemake.output.output_bam

iBAM = pysam.AlignmentFile(input_bam, 'rb')
oBAM = pysam.AlignmentFile(output_bam, 'wb', template=iBAM)

# Only flush winning pairs at chromosome boundaries if the input is coordinate
# sorted; otherwise duplicates of a family are not guaranteed to be contiguous.
sort_order = iBAM.header.to_dict().get('HD', {}).get('SO')
coord_sorted = (sort_order == 'coordinate')


def mean_qual(read):
    """Mean Phred base quality over the whole read, or None if unavailable."""
    quals = read.query_qualities
    if quals is None or len(quals) == 0:
        return None
    return sum(quals) / len(quals)


def flush_MB(MB, oBAM):
    for _, (_, read1, read2) in MB.items():
        oBAM.write(read1)
        oBAM.write(read2)


# Best pair seen so far per duplicate-family key: key -> (score, read1, read2)
MB = dict()
# Buffer for the first-seen mate of a pair, keyed by query_name
read_pairs = {}
current_chrom = None

for read in iBAM:
    # Skip secondary/supplementary alignments so a third record with the same
    # query_name can never be mis-paired with a primary mate.
    if read.is_secondary or read.is_supplementary:
        continue

    # Skip unpaired, unmapped, or unmated reads
    if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
        continue

    # Coordinate-sorted input: once we advance to a new chromosome, no further
    # pair can complete on the previous one, so flush its winners. Do NOT touch
    # read_pairs, since a buffered mate may still live on a later chromosome.
    if coord_sorted and read.reference_id != current_chrom:
        if current_chrom is not None:
            flush_MB(MB, oBAM)
            MB.clear()
        current_chrom = read.reference_id

    qname = read.query_name

    # Buffer until both mates are available
    if qname not in read_pairs:
        read_pairs[qname] = read
        continue

    mate = read_pairs.pop(qname)

    # Order the reads consistently
    read1, read2 = (read, mate) if not read.is_read2 else (mate, read)

    # Extract molecular barcode
    try:
        mb = read1.query_name.split('_MolecularBarcode:')[1]
    except IndexError:
        continue

    # Reference names
    chrom1 = iBAM.get_reference_name(read1.reference_id)
    chrom2 = iBAM.get_reference_name(read2.reference_id)

    # Strand-aware 5' start positions
    start1 = read1.reference_start if not read1.is_reverse else read1.reference_end
    start2 = read2.reference_start if not read2.is_reverse else read2.reference_end

    # Duplicate-family key uses info from both mates
    key = f"{chrom1}_{start1}_{read1.cigarstring}_{chrom2}_{start2}_{read2.cigarstring}_{mb}"

    # Combined mean base quality of the pair (skip pairs missing qualities)
    q1 = mean_qual(read1)
    q2 = mean_qual(read2)
    if q1 is None or q2 is None:
        continue
    score = q1 + q2

    # Keep the pair with the highest combined mean base quality
    if key not in MB or score > MB[key][0]:
        MB[key] = (score, read1, read2)

# Write both mates of every remaining winning pair
flush_MB(MB, oBAM)

iBAM.close()
oBAM.close()
