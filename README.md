# fastq2EZbakR

**fastq2EZbakR is under active development and should not be considered a stable distribution yet**

Reimplementation of the bam2bakR Snakemake workflow that makes several important improvements:

1. Better and more configurable preprocessing and alignment of fastq files
2. Fastq files are QCed with fastQC
3. Config file is more well organized and better documented
4. Automated unit testing has been established
5. Use of featureCounts for faster and more flexible read assignment to transcriptomic features (bam2bakR currently uses the slower HTSeq, which also doesn't properly account for splice junction mapping reads).
6. Lots of new functionality to support ongoing/soon-to-be-released work.

This pipeline will also be the site of exciting future development to expand the functionality of this pipeline, so stay tuned!

Documentation can be found here: https://fastq2ezbakr.readthedocs.io/en/latest/
