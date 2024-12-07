### GENERAL IDEA
# The input csv file has been sorted by the columns I would like to summarise by.
# All I need to do is count the number of rows with identical values for all
# of the sorted columns. The output row is then the values for the columns being kept
# as well as a new column, n, for the number of identical values observed.

import csv

import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    merged_table = snakemake.input
    cols_to_sum = snakemake.params.get("cols_to_sum")
    output_table = snakemake.output
    sample = snakemake.wildcards.sample
    
    ### Parse cols_to_sum

    cols_to_sum = cols_to_sum.split(',')

    header = ['sample'] + cols_to_sum + ['n']


    ### Create reader object for the merged table
    mt_f = open(merged_table[0], 'r')
    mt_r = csv.reader(mt_f)

    # Iterate past header
    full_header = next(mt_r)

    # Find which elements I want to summarise over
    sum_elements = [index for index, item in enumerate(full_header) if item in cols_to_sum]

    ### Create writer object for the output
    out_f = open(output_table[0], 'w', newline = '')
    out_w = csv.writer(out_f)
    out_w.writerow(header)

    ### Iterate over table
    first_row = True
    for row in mt_r:
        
        query = [row[i] for i in sum_elements]

        if first_row:

            # Initialize subject row to search for
            first_row = False
            subject = query
            count = 1
        
        elif query != subject:
            
            # Write current subject row + number of instances of that row found
            next_row = [sample] + subject + [count]
            out_w.writerow(next_row)
            count = 1

            # New subject row to search for
            subject = [row[i] for i in sum_elements]

        else:

            # Another identical row
            count += 1

    # Write the final row
    next_row = [sample] + subject + [count]
    out_w.writerow(next_row)

    ### Close connections
    out_f.close()
    mt_f.close()




