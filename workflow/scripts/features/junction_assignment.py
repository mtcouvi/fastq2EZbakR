'''
    File name: count_reads_only.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    Orcid: 
    Date created: 24/05/2022
    Date last modified: 24/05/2022
    Version: 1.0.0
    License: GPLv3
    Python Version: 
    Packages: pysam 
    Description: Script for creating table of junctions read mapped to
'''

import os
import pysam
import csv
import argparse
import datetime
import sys

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is a script to assign reads to junctions they overlap')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-b', '--bam', type=str, required=True, metavar = 'in_file.bam',
                    help='Bam file to process')
parser.add_argument('--filterBySJout', action='store_true',
                    help="Do not output data for reads overlapping junctions not included in STAR's SJout.tab file")
parser.add_argument('--filterNonCanon', action='store_true',
                    help="Do not output data for reads overlapping junctions that are non-canonical")
args = parser.parse_args()


# Get input bam file
inputName = args.bam.split('.bam')[0] 
output = inputName + '_junctions.csv'


# Assess filtering conditoin
SJfilter = args.filterBySJout
NCfilter = args.filterNonCanon

#  Set .csv file for writing (simulating _counts.rds file)
header = ['qname', 'junction_start', 'junction_end']
myfile = open(output, 'w', newline='')
wr = csv.writer(myfile)
wr.writerow(header)


# Row
r_info = [''] + 2*[0]

# Set .bam file for reading
samfile = pysam.AlignmentFile(args.bam, 'rb')



for read in samfile:


    junction_tags = [(tag, value) for tag, value in read.tags if tag in ['jI', 'jM']]

    if junction_tags[0][0] == 'jM':
        
        jM = junction_tags[0][1]
        jI = junction_tags[1][1]
    
    else:
        
        jM = junction_tags[1][1]
        jI = junction_tags[0][1]

    check = jM[0] != -1
    
    SJcheck = True
    NCcheck = True

    if SJfilter:
        SJcheck = jM[0] >= 20

    if NCfilter:
        NCcheck = jM[0] != 0 and jM[0] != 20

    if check and SJcheck and NCcheck:

        r_info[0] = read.query_name

        nj = len(jM)

        for j in range(nj):

            index = j*2


            if j == 0:

                r_info[1] = str(jI[index])
                r_info[2] = str(jI[index + 1])    

            else:

                r_info[1] = r_info[1] + "+" + str(jI[index])
                r_info[2] = r_info[2] + "+" + str(jI[index + 1])

            
        wr.writerow(r_info)
    

myfile.close()