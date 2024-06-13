### STRATEGY
# Step in pipeline whose RAM usage is a function of
# file size is the merge_features_and_muts rule. Also,
# cB files can similarly be created without using excessive
# RAM. To address these problems, I need to:
# 1) Sort the files to be merged by the read name column
# 2) Iterate through each file row by row and perform
# a left join with the mutation counts.
# 3) Sort files by columns that will be kept in cB
# 4) Iterate again this time counting the number of identical
# instances of the set of columns to keep. 

### Sorting rules

# Sort by read name (qname)
rule sort_mutcounts_by_qname:
    input:
        "results/counts/{sample}_counts.csv.gz",
    output:
        temp("results/sort_mutcounts_by_qname/{sample}_counts.csv.gz")
    shell:
        """
        sort -k1 -V {input} > {output}
        """

rule sort_fcgene_by_qname:
    input:
        "results/featurecounts_genes/{sample}.featureCounts",
    output:
        temp("results/sort_fcgene_by_qname/{sample}.featureCounts")
    shell:
        """
        sort -k1 -V {input} > {output}
        """


rule sort_fcexon_by_qname:
    input:
        "results/featurecounts_exons/{sample}.featureCounts",
    output:
        temp("results/sort_fcexon_by_qname/{sample}.featureCounts")
    shell:
        """
        sort -k1 -V {input} > {output}
        """


rule sort_junction_by_qname:


rule sort_fcee_by_qname:
    input:
        "results/featurecounts_ee/{sample}.featureCounts",
    output:
        temp("results/sort_fcee_by_qname/{sample}.featureCounts")
    shell:
        """
        sort -k1 -V {input} > {output}
        """


rule sort_fcei_by_qname:
    input:
        "results/featurecounts_ei/{sample}.featureCounts",
    output:
        temp("results/sort_fcei_by_qname/{sample}.featureCounts")
    shell:
        """
        sort -k1 -V {input} > {output}
        """