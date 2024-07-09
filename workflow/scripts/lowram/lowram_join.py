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


### SOME IDEAS
# Can pass as parameters booleans of whether or not a given feature assignment strategy was used
#
# Then can define the set of files to open and their respective names
# 
# Create csv/tsv reader object only for those files that exist
# 
# Loop over mutation counts and see which feature files have the qname and have assigned the read
# to an annotated feature.
#
# If present, new row will be all of the rows of mutation counts + the feature assignment. 
# NOTE: here I can even solve the featureCounts target list ordering problem myself!
#
# If absent, relevant feature assignment column will be "__no_feature"
#
# Write new row to output and loop on to the next

import csv

import sys

# String comparison with < yields incorrect results; this works
from natsort import natsort_keygen, ns
version_key = natsort_keygen(alg=ns.REAL)

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    input_files = snakemake.input
    output_table = snakemake.output


    ### See which file types are present

    genes_present = snakemake.params.get("genes")
    exons_present = snakemake.params.get("exons")
    junctions_present = snakemake.params.get("junctions")
    eb_present = snakemake.params.get("exonic_bins")
    ee_present = snakemake.params.get("ee")
    ei_present = snakemake.params.get("ei")
    transcripts_present = snakemake.params.get("transcripts")
    bft_present = snakemake.params.get("bamfile_transcripts")


    ### Open connections to files
    # Files will appear in this order:
    # 1) Mutation counts
    # 2) Genes
    # 3) Exons
    # 4) transcripts (feature counts)
    # 5) Exonic bins
    # 6) Bamfile transcripts
    # 7) Junctions
    # 8) eej
    # 9) eij

    mutcount = input_files[0]
    mutf = open(mutcount, 'r')
    mutr = csv.reader(mutf)
    count = 1

    header = next(mutr)

    if genes_present:
        genef = open(input_files[count], 'r')
        gener = csv.reader(genef, delimiter = '\t')
        count += 1
        header = header + ['GF']
        gene_row = next(gener)

    if exons_present:
        exonf = open(input_files[count], 'r')
        exonr = csv.reader(exonf, delimiter = '\t')
        count += 1
        header = header + ['XF']
        exon_row = next(exonr)


    if transcripts_present:
        transcriptf = open(input_files[count], 'r')
        transcriptr = csv.reader(transcriptf, delimiter = '\t')
        count += 1
        header = header + ['transcripts']
        transcript_row = next(transcriptr)


    if eb_present:
        ebf = open(input_files[count], 'r')
        ebr = csv.reader(ebf, delimiter = '\t')
        count += 1
        header = header + ['exon_bin']
        eb_row = next(ebr)


    if bft_present:
        bftf = open(input_files[count], 'r')
        bftr = csv.reader(bftf)
        count += 1
        header = header + ['bamfile_transcripts']

        # Iterate past header
        bft_row = next(bftr)
        bft_row = next(bftr)


    if junctions_present:
        jf = open(input_files[count], 'r')
        jr = csv.reader(jf)
        count += 1
        header = header + ['junction_start', 'junction_end']

        # Iterate past header
        j_row = next(jr)
        j_row = next(jr)


    if ee_present:
        eef = open(input_files[count], 'r')
        eer = csv.reader(eef, delimiter = '\t')
        count += 1
        header = header + ['ee_junction_id']
        ee_row = next(eer)


    if ei_present:
        eif = open(input_files[count], 'r')
        eir = csv.reader(eif, delimiter = '\t')
        count += 1
        header = header + ['ei_junction_id']
        ei_row = next(eir)

    ### Functions for handling the different scenarios

    def handle_featurecounts(query_row, iterator, current_outrow, mutrow):

        while query_row and version_key(query_row[0]) < version_key(mutrow[0]):
            query_row = next(iterator, None)

        if query_row and query_row[0] == mutrow[0]:

            if query_row[1] == "Assigned":

                additional_out = query_row[3]

                if ',' in additional_out:

                    additional_out = additional_out.split(',')
                    additional_out = sorted(additional_out)
                    additional_out = '+'.join(additional_out)

                current_outrow = current_outrow + [additional_out]

            else:

                current_outrow = current_outrow + ['__no_feature']    

        else:

            current_outrow = current_outrow + ['__no_feature']

        return current_outrow, query_row, iterator


    def handle_bamfilet(query_row, iterator, current_outrow, mutrow):

        while query_row and version_key(query_row[0]) < version_key(mutrow[0]):
            query_row = next(iterator, None)

        if query_row and query_row[0] == mutrow[0]:

            current_outrow = current_outrow + [query_row[1]]

        else:

            current_outrow = current_outrow + ['__no_feature']

        return current_outrow, query_row, iterator


    def handle_junctions(query_row, iterator, current_outrow, mutrow):

        while query_row and version_key(query_row[0]) < version_key(mutrow[0]):
            query_row = next(iterator, None)

        if query_row and query_row[0] == mutrow[0]:

            current_outrow = current_outrow + [query_row[1]] + [query_row[2]]

        else:

            current_outrow = current_outrow + ['__no_feature', '__no_feature']

        return current_outrow, query_row, iterator


    ### Create merged output

    # Open sorted files and output file
    with open(output_table[0], 'w', newline='') as output_file:

        writer = csv.writer(output_file)
        
        # Get the header and write to the output
        writer.writerow(header)
        
        # Initialize variables
        for row_m in mutr:

            outrow = row_m
                
            # Add GF information
            if genes_present:

                outrow, gene_row, gener = handle_featurecounts(gene_row, gener, outrow, row_m)

            # Add XF information
            if exons_present:

                outrow, exon_row, exonr = handle_featurecounts(exon_row, exonr, outrow, row_m)

            # Add transcripts information
            if transcripts_present:

                outrow, transcript_row, transcriptr = handle_featurecounts(transcript_row, transcriptr, outrow, row_m)

            # Add transcripts information
            if eb_present:

                outrow, eb_row, ebr = handle_featurecounts(eb_row, ebr, outrow, row_m)

            # Add transcripts information
            if bft_present:

                outrow, bft_row, bftr = handle_bamfilet(bft_row, bftr, outrow, row_m)
            
            # Add transcripts information
            if junctions_present:

                outrow, j_row, jr = handle_junctions(j_row, jr, outrow, row_m)
            
            if ee_present:

                outrow, ee_row, eer = handle_featurecounts(ee_row, eer, outrow, row_m)
            
            if ei_present:

                outrow, ei_row, eir = handle_featurecounts(ei_row, eir, outrow, row_m)
            
            writer.writerow(outrow)


    ### Close file connections

    mutf.close()

    if genes_present:
        genef.close()

    if exons_present:
        exonf.close()

    if transcripts_present:
        transcriptf.close()

    if eb_present:
        ebf.close()

    if bft_present:
        bftf.close()

    if junctions_present:
        jf.close()

    if ee_present:
        eef.close()

    if ei_present:
        eif.close()
