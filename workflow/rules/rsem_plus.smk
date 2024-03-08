# Convert RSEM bam file to a csv with transcript probabilities for RSEM+
rule rsem_to_csv:
    input:
        "results/rsem/{sample}.transcript.bam",
    output:
        "results/rsem_csv/{sample}_rsem.csv.gz",
        temp("results/rsem_csv/{sample}_check.txt"),
    params:
        shellscript=workflow.source_path("../scripts/rsem_plus/rsem_to_csv.sh"),
        pythonscript=workflow.source_path("../scripts/rsem_plus/rsem_csv.py"),
        awkscript=workflow.source_path("../scripts/rsem_plus/fragment_sam_rsem.awk"),
    log:
        "logs/rsem_to_csv/{sample}.log",
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.pythonscript} {params.awkscript} 1> {log} 2>&1
        """


# Estimate transcript isoform fraction news
rule transcript_fn:
    input:
        rsem="results/rsem_csv/{sample}_rsem.csv.gz",
        counts="results/merge_features_and_muts/{sample}_counts.csv.gz",
    output:
        outfile="results/transcript_fn/{sample}_RSEM_plus.csv",
    params:
        rscript=workflow.source_path("../scripts/rsem_plus/RSEM_plus.R"),
        pnew=get_pnew,
        pold=get_pold,
    log:
        "logs/transcript_fn/{sample}.log",
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        r"""
        chmod +x {params.rscript}
        {params.rscript} -o {output.outfile} -c {input.counts} -r {input.rsem} -s {wildcards.sample} -n {params.pnew} -b {params.pold} 1> {log} 2>&1
        """


# Combine transcript isoform fraction new estimates for all samples
## THIS IS WRONG AND PUTS HEADERS INSIDE THE MIDDLE OF THE TABLE
## Thankfully read_csv() is able to properly parse this shit
rule combine_fn:
    input:
        expand("results/transcript_fn/{samps}_RSEM_plus.csv", samps=SAMP_NAMES),
    output:
        "results/transcript_fn/RSEM_plus.csv",
    log:
        "logs/combine_fn/combine_fn.log",
    threads: 1
    conda:
        "../envs/full.yaml"
    shell:
        """
        cat {input} >> {output}
        """
