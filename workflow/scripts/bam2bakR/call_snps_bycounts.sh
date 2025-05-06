#!/bin/bash
set -euo pipefail

# Cute trick to deal with fact that there is an uncertain number of control sample
    # control_samples becomes array with all args
    # remove args that I know aren't the actual control_samples

cpus=$1
nsamps=$2
output_txt=$3
snp_counts=$4
genome_fasta=$5

shift 8

control_samples=("$@")




if [[ "$nsamps" -gt 0 ]] 
then
    # Loop through control samples:
        declare -a bam_list=()

        for cs in ${control_samples[@]}; do
            name=$(echo "$cs" | cut -d '/' -f 3 | rev | cut -c8- | rev)

            sorted="./results/snps/"$name"_sort.bam"

            samtools sort -@ "$cpus" -o ./results/snps/"$name"_sort.bam "$cs"
            samtools index -@ "$cpus" ./results/snps/"$name"_sort.bam
            #bcftools mpileup --threads "$cpus" -f "$genome_fasta" ./results/snps/"$name"_sort.bam | bcftools call --threads "$cpus" -mv > ./results/snps/snp-"$name".vcf

            bam_list+=( "$sorted" )

        done


        echo "${bam_list[@]}" | tr ' ' '\n' > ./results/snps/bam.list

        bcftools mpileup --threads "$cpus" \
                         -f "$genome_fasta" \
                         -b ./results/snps/bam.list \
                         -a AD,DP \
                         -Ou \
        | bcftools view --threads "$cpus" -i 'FORMAT/AD[0:1]>=${snp_counts}' -o Min${snp_counts}_sites.vcf


		# Make snp.txt
		grep -v 'INDEL' Min${snp_counts}_sites.vcf > Min${snp_counts}_sites_noINDEL.vcf
		awk '{if($1 !~ /^#/){split($5,alt,","); print $4 ":" alt[1] ":" $1 ":" $2}}' Min${snp_counts}_sites_noINDEL.vcf | sort | uniq > $output_txt
        
        echo '* SNPs called and snp.txt generated'

else
    touch "$output_txt"
fi


rm -f 0