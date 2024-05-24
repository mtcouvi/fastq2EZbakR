### THESE RULES PERTAIN TO THE ASSIGNMENT OF READS TO FEATURES WITH FEATURECOUNTS


# Assign reads to genes
rule featurecounts_genes:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_genes/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_genes_extra"] + FC_GENES_PARAMS,
    log:
        "logs/featurecounts_gene/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Assign reads to exons
rule featurecounts_exons:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_exons/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_exons_extra"] + FC_EXONS_PARAMS,
    log:
        "logs/featurecounts_exons/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Assign reads to transcripts
rule featurecounts_transcripts:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_transcripts/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_transcripts_extra"] + FC_TRANSCRIPTS_PARAMS,
    log:
        "logs/featurecounts_transcripts/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Assign reads to exonic bins
rule featurecounts_exonbins:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["flat_annotation"],
    output:
        multiext(
            "results/featurecounts_exonbins/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_exonbins_extra"] + FC_EXONBINS_PARAMS,
    log:
        "logs/featurecounts_exonbins/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Get the set of isoforms a read maps to from the transcriptome bam
# TO-DO: No reason this can't be split up and multi-threaded
rule read_to_transcripts:
    input:
        bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
    output:
        table="results/read_to_transcripts/{sample}.csv",
    log:
        "logs/read_to_transcripts/{sample}.log",
    conda:
        "../envs/full.yaml"
    threads: 1
    script:
        "../scripts/features/transcript_assignment.py"


# Get set of junctions a read overlaps
rule read_to_junctions:
    input:
        "results/align/{sample}.bam",
    output:
        "results/read_to_junctions/{sample}.csv.gz",
        temp("results/read_to_junctions/{sample}_check.txt"),
    params:
        shellscript=workflow.source_path("../scripts/features/junction_assignment.sh"),
        pythonscript=workflow.source_path("../scripts/features/junction_assignment.py"),
        awkscript=workflow.source_path("../scripts/rsem_plus/fragment_sam_rsem.awk"),
    log:
        "logs/read_to_junctions/{sample}.log",
    threads: 32
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.pythonscript} {params.awkscript} 1> {log} 2>&1
        """
