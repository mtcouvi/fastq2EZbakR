__author__ = "Isaac Vock"
__copyright__ = "Copyright 2023, Isaac Vock"
__email__ = "isaac.vock@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
import os.path as path
import sys

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
n = len(snakemake.input.fastqs)

assert(
    n == 1 or n == 2
), "input->fastqs must have 1 (single-end) or 2 (paired-end) elements"

fastqs = snakemake.input.fastqs
output = snakemake.output.unzipped_fqs

if n == 1:

    shell(
        "(pigz -d -p {snakemake.threads} -c {fastqs} > {output}) {log}"
    )   

else:

    fastq1 = fastqs[0]
    fastq2 = fastqs[1]
    out1 = output[0]
    out2 = output[1]

    shell(
        "(pigz -d -p {snakemake.threads} -c {fastq1} > {out1}) {log}"
    ) 

    shell(
        "(pigz -d -p {snakemake.threads} -c {fastq2} > {out2}) {log}"
    )  


