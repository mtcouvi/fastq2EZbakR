#!/bin/bash

# Source the paths and variables:
cpus=$1
pos_cutoff=$2
mutposout=$3
mutposfilter=$4
high_cutoff=$5



# Read all _cU.csv.gz files and save them as cU-DATE.csv.gz
parallel -j $cpus --plus "cat <(echo Filename:{1%_cU.csv.gz}) <(pigz -d -k -c -p 1 {1})" ::: ./results/counts/*_cU.csv.gz \
    | awk -v OFS="," '
            $1 ~ /Filename/ {
                split($1, sample, ":")
                next
            }
            NR == 2 {
                header = $0
                print "sample", $0
                next
            }
            $0 == header { 
                next
            }
            {
                print sample[2], $0
            }' \
    | awk '{ if (NR > 1) {$1 = substr($1, 18); print } else print }' \
    | pigz -p $cpus > "$mutposout"


pigz -d -c "$mutposout" \
| awk -F "," \
        -v cutoff="$pos_cutoff" \
        -v upper="$high_cutoff" \
        'NR == 1 || $(NF-1) >= cutoff && $(NF-1) <= upper { print}' | pigz -p $cpus >  "$mutposfilter"


echo "**  site-specific mutation file created: mutpos.csv.gz"


