# fastq2EZbakR

**fastq2EZbakR is under active development and should not be considered a stable distribution yet**

Reimplementation of the bam2bakR Snakemake workflow that makes several important improvements:

1. Better and more configurable preprocessing and alignment of fastq files
2. Fastq files are QCed with fastQC
3. Config file is more well organized and better documented
4. Automated unit testing has been established

This pipeline will also be the site of exciting future development to expand the functionality of this pipeline, so stay tuned!

The steps for running this pipline are identical to [bam2bakR](https://tl-snakemake.readthedocs.io/en/latest/), though there are differences in the config file to be aware of. In addition, if you are a Yale user or are interested in seeing how to optimize deployment of the pipeline on a shared cluster (in particular one which uses the slurm workload manager), check out the documentation for how to run a separate Snakemake pipeline I developed called [PROseq_etal](https://proseq-etal.readthedocs.io/en/latest/simon/). The steps for deploying PROseq_etal are identical to that of deploying fastq2EZbakR, and the same profile can be used to manage job scheduling in both cases.
