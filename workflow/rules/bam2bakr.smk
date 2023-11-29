### BIG PICTURE TO-DO
## 1) Make compatible with new config parameter names
## 2) Change names of rules and output directories to make more sense

### TO-DO
## 1) Clean up sort/filter function; maybe reduce number of output files
if config["bam2bakr"]:

    if config["remove_tags"]:
        
        # Remove tags from bam files that can break HTSeq
        rule remove_tags:
            input:
                input_bam=get_input_bams,
            output:
                output_bam="results/remove_tags/{sample}_no_jI_jM.bam",
            log:
                "logs/remove_tags/{sample}.log"
            conda:
                "../envs/full.yaml"
            script:
                "../scripts/bam2bakR/remove_tags.py"

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                "results/remove_tags/{sample}_no_jI_jM.bam",
            output:
                "results/sf_reads/{sample}.s.sam",
                "results/sf_reads/{sample}_fixed_mate.bam",
                "results/sf_reads/{sample}.f.sam",
            log:
                "logs/sort_filter/{sample}.log"
            params: 
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
                format=FORMAT
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1
                """

    else:

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                get_input_bams
            output:
                "results/sf_reads/{sample}.s.sam",
                "results/sf_reads/{sample}_fixed_mate.bam",
                "results/sf_reads/{sample}.f.sam",
            log:
                "logs/sort_filter/{sample}.log"
            params: 
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
                format=FORMAT
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1
                """


else:

    # Filter out multi-mappers and sort reads
    rule sort_filter:
        input:
            "results/align/{sample}.bam"
        output:
            "results/sf_reads/{sample}.s.sam",
            "results/sf_reads/{sample}_fixed_mate.bam",
            "results/sf_reads/{sample}.f.sam",
        log:
            "logs/sort_filter/{sample}.log"
        params: 
            shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
            format=FORMAT
        threads: 8
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1
            """



### TO-DO:
## 1) Chunking and parallel processing of bam files
## 2) Allow users to specify various parameters
# Use custom htseq script to quantify features 
# Also creates bam files with tag designating feature that each read was mapped to; useful during mutation counting
if config["strategies"]["FlatStacks"]:

    rule htseq_cnt:
        input:
            sam="results/sf_reads/{sample}.s.sam",
            flatstack=config["flat_annotation"]
        output:
            "results/htseq/{sample}_tl.bam",
            temp("results/htseq/{sample}_check.txt")
        params: 
            shellscript=workflow.source_path("../scripts/bam2bakR/htseq.sh"),
            pythonscript=workflow.source_path("../scripts/bam2bakR/count_triple.py"),
            strand=config["strandedness"],
            flattened=config["strategies"]["FlatStacks"],
        log:
            "logs/htseq_cnt/{sample}.log"
        threads: 3
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            chmod +x {params.pythonscript}
            {params.shellscript} {threads} {wildcards.sample} {input.sam} {output} {input.flatstack} {params.strand} {params.pythonscript} {params.flattened} 1> {log} 2>&1
            """

else:

    rule htseq_cnt:
        input:
            sam="results/sf_reads/{sample}.s.sam",
            annotation=config["annotation"]
        output:
            "results/htseq/{sample}_tl.bam",
            temp("results/htseq/{sample}_check.txt")
        params: 
            shellscript=workflow.source_path("../scripts/bam2bakR/htseq.sh"),
            pythonscript=workflow.source_path("../scripts/bam2bakR/count_triple.py"),
            strand=config["strandedness"],
            flattened=config["strategies"]["FlatStacks"],
        log:
            "logs/htseq_cnt/{sample}.log"
        threads: 3
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            chmod +x {params.pythonscript}
            {params.shellscript} {threads} {wildcards.sample} {input.sam} {output} {input.annotation} {params.strand} {params.pythonscript} {params.flattened} 1> {log} 2>&1
            """

### TO-DO
## 1) Properly log standard out
# Calculate normalization scale factor to be applied to tracks        
if NORMALIZE:
    rule normalize:
        input:
            expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
        output:
            "results/normalization/scale"
        log:
            "logs/normalize/normalize.log"
        params:
            rscript=workflow.source_path("../scripts/bam2bakR/normalize.R"),
            spikename=config["spikename"]
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            r"""
            chmod +x {params.rscript}
            {params.rscript} --dirs ./results/htseq/ --spikename {params.spikename}
            mv scale {output}
            """
else:
    rule normalize:
        input:
            expand("results/htseq/{sample}_tl.bam", sample = SAMP_NAMES)
        output:
            "results/normalization/scale"
        log:
            "logs/normalize/normalize.log"
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            touch {output}
            """

# Index genome fasta file for snp calling
rule genome_index:
    input:
        str(config["genome"])
    output:
        get_index_name()
    log:
        "logs/genome_index/genome-faidx.log",
    threads: 1
    conda:
        "../envs/index.yaml"
    script:
        "../scripts/bam2bakR/genome-faidx.py"

## TO-DO
# 1) Allow users to provide custom SNP file
# Identify SNPs to be accounted for when counting mutations
rule call_snps:
    input:
        str(config["genome"]),
        get_index_name(),
        expand("results/htseq/{ctl}_tl.bam", ctl = CTL_NAMES)
    params:
        nctl = nctl,
        shellscript = workflow.source_path("../scripts/bam2bakR/call_snps.sh"),
    output:
        "results/snps/snp.txt",
        "results/snps/snp.vcf",
        temp("results/snps/mkdir.txt")
    log:
        "logs/call_snps/ctl_samps.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {params.nctl} {output} {input} 1> {log} 2>&1
        """

# TO-DO:
# 1) Add mutation position optimizations and functionality
# Count mutations 
rule cnt_muts:
    input:
        "results/htseq/{sample}_tl.bam",
        "results/snps/snp.txt"
    params:
        format = FORMAT,
        minqual = config["minqual"],
        mut_tracks = config["mut_tracks"],
        strand = STRAND,
        shellscript = workflow.source_path("../scripts/bam2bakR/mut_call.sh"),
        pythonscript = workflow.source_path("../scripts/bam2bakR/mut_call.py"),
        awkscript = workflow.source_path("../scripts/bam2bakR/fragment_sam.awk")
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.minqual} {params.mut_tracks} {params.format} {params.strand} {params.pythonscript} {params.awkscript} 1> {log} 2>&1
        """

### Get the set of transcripts that a read aligned to
### and combine that info with mutation counting on genome
### aligned bams
if config["strategies"]["Transcripts"]:

    rule read_to_transcripts:
        input:
            bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
        output:
            table=temp("results/read_to_transcripts/{sample}.csv")
        log:
            "logs/read_to_transcripts/{sample}.log"
        conda:
            "../envs/full.yaml"
        threads: 1
        script:
            "../scripts/bam2bakR/count_transcriptome.py"

    # Sort transcript mapping table
    # to facilitate memory efficient joining later
    rule sort_transcripts_table:
        input:
            transcripts="results/read_to_transcripts/{sample}.csv",
        output:
            sorted="results/read_to_transcripts/{sample}_sorted.csv"
        log:
            "logs/sort_transcripts_table/{sample}.log"
        params:
            script = workflow.source_path("../scripts/bam2bakR/cheap_sort.sh"),
            lines = CHUNK_SIZE,
        conda:
            "../envs/full.yaml"
        threads: 1
        shell:
            """
            chmod +x {params.script}
            {params.script} {input} {output} 1 {params.lines} FALSE ./results/read_to_transcripts {wildcards.sample}
            """
    
    # Sort mutation counts
    # to facilitate memory efficient joining later
    rule sort_counts:
        input:
            counts="results/counts/{sample}_counts.csv.gz",
        output:
            sorted="results/sort_counts/{sample}_sorted.csv"
        log:
            "logs/sort_counts/{sample}.log"
        params:
            script = workflow.source_path("../scripts/bam2bakR/cheap_sort.sh"),
            lines = CHUNK_SIZE,
        conda:
            "../envs/full.yaml"
        threads: 1
        shell:
            """
            chmod +x {params.script}
            {params.script} {input} {output} 1 {params.lines} TRUE ./results/sort_counts {wildcards.sample}
            """

    rule merge_counts:
        input:
            muts="results/sort_counts/{sample}_sorted.csv",
            transcripts="results/read_to_transcripts/{sample}_sorted.csv"
        output:
            merged=temp("results/merge_counts/{sample}_counts.csv")
        conda:
            "../envs/full.yaml"
        threads: 1
        script:
            "../scripts/bam2bakR/noram_join.py"

    rule compress_counts:
        input:
            "results/merge_counts/{sample}_counts.csv"
        output:
            "results/merge_counts/{sample}_counts.csv.gz"
        conda:
            "../envs/full.yaml"
        threads: 8
        shell:
            "pigz -p {threads} -c {input} > {output}"


    # Make cB file that will be input to bakR
    rule makecB:
        input:
            expand("results/merge_counts/{sample}_counts.csv.gz", sample=SAMP_NAMES)
        output:
            "results/cB/cB.csv.gz"
        params:
            shellscript = workflow.source_path("../scripts/bam2bakR/master.sh"),
            keepcols = config["keepcols"],
            mut_tracks=config["mut_tracks"]
        log:
            "logs/makecB/master.log"
        threads: 20
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {output} {params.keepcols} {params.mut_tracks} ./results/merge_counts TRUE 1> {log} 2>&1
            """


else:

    # Make cB file that will be input to bakR
    rule makecB:
        input:
            expand("results/counts/{sample}_counts.csv.gz", sample=SAMP_NAMES)
        output:
            "results/cB/cB.csv.gz"
        params:
            shellscript = workflow.source_path("../scripts/bam2bakR/master.sh"),
            keepcols=config["keepcols"],
            mut_tracks=config["mut_tracks"]
        log:
            "logs/makecB/master.log"
        threads: 20
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {output} {params.keepcols} {params.mut_tracks} ./results/counts FALSE 1> {log} 2>&1
            """


# Make color-coded tracks
rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/htseq/{sample}_tl.bam",
	    "results/normalization/scale"
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=config["mut_tracks"], id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    params:
        shellscript = workflow.source_path("../scripts/bam2bakR/tracks.sh"),
        pythonscript = workflow.source_path("../scripts/bam2bakR/count_to_tracks.py"),
        mut_tracks=config["mut_tracks"],
        genome=config["genome"],
        WSL=config["WSL"],
        normalize=config["normalize"]
    log:
        "logs/maketdf/{sample}.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {params.mut_tracks} {params.genome} {params.WSL} {params.normalize} {params.pythonscript} {output} 1> {log} 2>&1
        """

