# fastq2EZbakR

**fastq2EZbakR is under active development and should not be considered a stable distribution yet**

Reimplementation of the bam2bakR Snakemake workflow that makes several important improvements:

1. Better and more configurable preprocessing and alignment of fastq files
2. Fastq files are QCed with fastQC
3. Config file is more well organized and better documented
4. Automated unit testing has been established
5. Use of featureCounts for faster and more flexible read assignment to transcriptomic features (bam2bakR currently uses the slower HTSeq, which also doesn't properly account for splice junction mapping reads).

This pipeline will also be the site of exciting future development to expand the functionality of this pipeline, so stay tuned!

The steps for running this pipline are identical to [bam2bakR](https://tl-snakemake.readthedocs.io/en/latest/), though there are differences in the config file to be aware of. In addition, if you are a Yale user or are interested in seeing how to optimize deployment of the pipeline on a shared cluster (in particular one which uses the slurm workload manager), check out the [documentation](https://pipelinedocs.readthedocs.io/en/latest/simon/) for how to generally run all of the Snakemake pipelines I have developed. The steps for deploying these pipelines are identical to that of deploying fastq2EZbakR, and the same profile can be used to manage job scheduling in both cases.
