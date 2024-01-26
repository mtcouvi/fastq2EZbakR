# Make transcriptome fasta from genome fasta + gtf
rule make_transcriptome_fasta:
    input:
        fasta=config["genome"],
        annotation=AandQ_ANNOTATION,
    output:
        records="make_transcriptome_fasta/transcriptome.fasta",
    threads: 1
    log:
        "logs/make_transcriptome_fasta/gffread.log",
    params:
        extra=config["gffread_extra"],
    conda:
        "../envs/gffread.yaml"
    script:
        "../scripts/gffread.py"


### Make Salmon indices

if config["decoy_settings"]["entire_genome"]:

    rule get_decoys:
        input:
            genome=config["genome"],
            transcriptome="make_transcriptome_fasta/transcriptome.fasta",
        output:
            gentrome="results/salmon_decoys/gentrome.fa",
            decoys="results/salmon_decoys/decoys.txt",
        log:
            "logs/get_decoys/salmon_decoys.log",
        threads: 2
        wrapper:
            "v2.6.0/bio/salmon/decoys"

else:

    rule get_decoys:
        input:
            genome=config["genome"],
            annotation=config["annotation"],
            transcriptome="make_transcriptome_fasta/transcriptome.fasta",
        output:
            gentrome="results/salmon_decoys/gentrome.fa",
            decoys="results/salmon_decoys/decoys.txt",
        log:
            "logs/get_decoys/salmon_decoys.log",
        threads: 8
        params:
            extra=config["salmon_generateDecoy_params"],
            shellscript=workflow.source_path("../scripts/generateDecoyTranscriptome.sh"),
        conda:
            "../envs/decoys.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {params.extra} -a {input.annotation} -o results/salmon_decoys/ \
            -j {threads} -g {input.genome} -t {input.transcriptome} 2> {log}
            """


if config["decoy_settings"]["make_decoy"]:

    rule index:
        input:
            sequences=SALMON_TRANSCRIPTOME,
            decoys=SALMON_DECOYS,
        output:
            multiext(
                config["salmon_indices"],
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
        log:
            "logs/index/salmon_index.log",
        threads: 16  # Borrowed from vignette here: https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/
        params:
            extra=config["salmon_index_params"],
        wrapper:
            "v2.6.0/bio/salmon/index"

else:

    rule index:
        input:
            sequences=SALMON_TRANSCRIPTOME,
        output:
            multiext(
                config["salmon_indices"],
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
        log:
            "logs/index/salmon_index.log",
        threads: 16  # Borrowed from vignette here: https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/
        params:
            extra=config["salmon_index_params"],
        wrapper:
            "v2.6.0/bio/salmon/index"


### Quantify with Salmon

if config["PE"]:

    rule quant:
        input:
            r1=get_fastq_r1,
            r2=get_fastq_r2,
            index=multiext(
                config["salmon_indices"],
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
        output:
            quant="results/quant/{sample}/quant.sf",
            lib="results/quant/{sample}/lib_format_counts.json",
        log:
            "logs/quant/{sample}_salmon.log",
        params:
            libtype=LIBTYPE,
            extra=config["salmon_quant_params"],
        threads: 12  # See https://salmon.readthedocs.io/en/latest/salmon.html Note for motivation
        wrapper:
            "v2.6.0/bio/salmon/quant"

else:

    rule quant:
        input:
            r=get_fastq_r1,
            index=multiext(
                config["salmon_indices"],
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
        output:
            quant="results/quant/{sample}/quant.sf",
            lib="results/quant/{sample}/lib_format_counts.json",
        log:
            "logs/quant/{sample}_salmon.log",
        params:
            libtype=LIBTYPE,
            extra=config["salmon_quant_params"],
        threads: 12
        wrapper:
            "v2.6.0/bio/salmon/quant"
