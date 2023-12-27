import glob
import os


### If RSEM+ is set, make sure that exons are quantified
if config["strategies"]["RSEMp"]:

    config["features"]["exons"] = True


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Sample names to help expanding lists of all bam files
# and to aid in defining wildcards
SAMP_NAMES = list(config['samples'].keys())


# Directory containing index; used in case of certain aligners
INDEX_DIR = config["indices"]


# Make life easier for users and catch if they add a '/' at the end of their path
# to alignment indices. If so, remove it to avoid double '/' 
if config["indices"].endswith('/'):
    INDEX_PATH = str(config["indices"])
    INDEX_PATH = INDEX_PATH[:-1]
else:
    INDEX_PATH = str(config["indices"])


# Determine how many fastqs to look for
if config["PE"]:
    READS = [1, 2]
    READ_NAMES = ['r1', 'r2']
else:
    READS = [1]
    READ_NAMES = ['r1']

# Get input fastq files for first step
def get_input_fastqs(wildcards):
    fastq_path = config["samples"][wildcards.sample]
    fastq_files = sorted(glob.glob(f"{fastq_path}/*.fastq*"))
    return fastq_files


# Check if fastq files are gzipped
fastq_paths = config["samples"]

is_gz = False

for p in fastq_paths.values():

    fastqs = sorted(glob.glob(f"{p}/*.fastq*"))
    test_gz = any(path.endswith('.fastq.gz') for path in fastqs)
    is_gz = any([is_gz, test_gz])



# Columns in final cB
keepcols = ['sample', 'sj', 'rname']

if config["features"]["genes"]:

    keepcols.append("GF")

if config["features"]["exons"]:

    keepcols.append("XF")

if config["features"]["transcripts"]:

    keepcols.append("transcripts")

if config["features"]["exonic_bins"]:

    keepcols.append("exon_bin")

keepcols = ','.join(keepcols)



### STAR HELPERS

# Trimmed fastq file paths, used as input for aligners

def get_fastq_r1(wildcards):
    if config["PE"]:

        return expand("results/trimmed/{SID}.1.fastq", SID = wildcards.sample)

    else:

        return expand("results/trimmed/{SID}.fastq", SID = wildcards.sample)

def get_fastq_r2(wildcards):
    if config["PE"]:

        return expand("results/trimmed/{SID}.2.fastq", SID = wildcards.sample)

    else:

        return ""


### Extra parameters passed to STAR alignment


STAR_PARAMS = str(config["star_align_params"])

## Add necessary tags and change other settings as necessary

args = STAR_PARAMS.split()

# Tags to always include
    # Technically, NH, HI, AS, and nM are added by default
    # NH = number of places the read aligns
    # HI = Numerical index that is useful to separate multi-mapping reads
    # NM = Edit distance between read and reference
    # AS = Alignment score
    # MD = Used for mutation counting
    # nM = Number of mismatches
sam_attributes = set(config["star_sam_tags"])

# Process --outSAMattributes
if "--outSAMattributes" in args:
    
    index = args.index("--outSAMattributes")

    # Assuming the attributes are space-separated and continuous after the flag
    next_flag_index = len(args)
    for i in range(index + 1, len(args)):
        if args[i].startswith('--'):
            next_flag_index = i
            break
    existing_attributes = set(args[index + 1:next_flag_index])
    combined_attributes = existing_attributes.union(sam_attributes)
    args[index + 1:next_flag_index] = list(combined_attributes)

else:
    args.extend(["--outSAMattributes"] + list(sam_attributes))

# Process --quantMode
quant_mode = set(["TranscriptomeSAM", "GeneCounts"])

if "--quantMode" in args:
    
    index = args.index("--quantMode")

    # Figure out what arguments were supplied to --quantMode
    next_flag_index = len(args)
    for i in range(index + 1, len(args)):
        if args[i].startswith('--'):
            next_flag_index = i
            break

    existing_attributes = set(args[index + 1:next_flag_index])

    # Add desired quantification modes
    combined_attributes = existing_attributes.union(quant_mode)
    args[index + 1:next_flag_index] = list(combined_attributes)

else:
    args.extend(["--quantMode"] + list(quant_mode))

# Force --outSAMtype to be BAM SortedByCoordinate
if "--outSAMtype" in args:
    
    index = args.index("--outSAMtype")

    # Replace the existing outSAMtype values with the desired ones
    args[index + 1:index + 3] = ["BAM", "SortedByCoordinate"]

else:
    
    args.extend(["--outSAMtype", "BAM", "SortedByCoordinate"])

# Force --sjdbGTFfile to be the provided annotation
if "--sjdbGTFfile" in args:
    
    index = args.index("--sjdbGTFfile")

    # Replace the existing outSAMtype values with the desired ones
    args[index + 1:index + 3] = [str(config["annotation"])]

else:
    
    args.extend(["--sjdbGTFfile", str(config["annotation"])])



STAR_EXTRA = ' '.join(args)


### HISAT2 HELPERS

# Figure out what to pass to --rna-strandedness
if config["strandedness"] == "no":

    HISAT2_STRANDEDNESS = ""

elif config["strandedness"] == "reverse":

    if config["PE"]:

        HISAT2_STRANDEDNESS = "--rna-strandness RF"

    else:

        HISAT2_STRANDEDNESS = "--rna-strandness R"



else:

    if config["PE"]:

        HISAT2_STRANDEDNESS = "--rna-strandness FR"

    else:

        HISAT2_STRANDEDNESS = "--rna-strandness F"

# Index base name
if config["hisat2_index_base"]:

    HISAT2_BASE = "{}/{}".format(INDEX_PATH, config["hisat2_index_base"])

else:

    HISAT2_BASE = "{}/{}".format(INDEX_PATH, "hisat2_index")





### BAM2EZBAKR PARAMETERS


# -s4U control sample names
CTL_NAMES = list(config['control_samples'])


# +s4U sample names
s4U_SAMPS = list(set(SAMP_NAMES) - set(CTL_NAMES))


# Number of -s4U control samples
nctl = len(CTL_NAMES)


# Number of +s4U samples
nsamps = len(SAMP_NAMES)


# Name for genome index
def get_index_name():
    genome = config["genome"]
    index = str(genome) + ".fai"
    return index


# PE or SE? 

if config["PE"]:

    FORMAT = "PE"

else:

    FORMAT = "SE"


# Calculate scale factors for tracks?
NORMALIZE = config['normalize']



# bam2EZbakR bam file fetching
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# pnew estimates for RSEM+
def get_pnew(wildcards):
    return config["pnews"][wildcards.sample]

# pold estimates for RSEM+
def get_pold(wildcards):
    return config["polds"][wildcards.sample]



### FEATURECOUNTS HELPERS
 
## All of the files to merge
def get_merge_input(wildcards):

    MERGE_INPUT = []

    MERGE_INPUT.append(expand("results/counts/{sid}_counts.csv.gz", sid = wildcards.sample))

    if config["features"]["genes"]:

        MERGE_INPUT.append(expand("results/featurecounts_genes/{sid}.featureCounts", sid = wildcards.sample))

    if config["features"]["exons"]:

        MERGE_INPUT.append(expand("results/featurecounts_exons/{sid}.featureCounts", sid = wildcards.sample))


    if config["features"]["transcripts"]:

        MERGE_INPUT.append(expand("results/featurecounts_transcripts/{sid}.featureCounts", sid = wildcards.sample))


    if config["features"]["exonic_bins"]:

        MERGE_INPUT.append(expand("results/featurecounts_exonbins/{sid}.featureCounts", sid = wildcards.sample))

    return MERGE_INPUT



# Get strandedness parameter
if config["strandedness"] == "reverse":

    FC_STRAND = 2

elif config["strandedness"] == "yes":

    FC_STRAND = 1

else:

    FC_STRAND = 0


## Get extra parameters for gene calling

if config["PE"]:

    FC_GENES_PARAMS = " -R -f -g gene_id -t transcript -p --countReadPairs"

else:

    FC_GENES_PARAMS = " -R -f -g gene_id -t transcript"


## Get extra parameters for transcript calling

if config["PE"]:

    FC_TRANSCRIPTS_PARAMS = " -R -f -g transcript_id -t exon -O -p --countReadPairs"

else:

    FC_TRANSCRIPTS_PARAMS= " -R -f -g transcript_id -t exon -O"


## Get extra parameters for exon bin calling


# Add paired-end status
if config["PE"]:

    FC_EXONBINS_PARAMS = " -R -f -g exon_id -t exonic_part -O -p --countReadPairs"

else:

    FC_EXONBINS_PARAMS= " -R -f -g exon_id -t exonic_part -O"



## Get extra parameters for exon calling


# Add paired-end status
if config["PE"]:

    FC_EXONS_PARAMS = " -R -g gene_id -J -p --countReadPairs"

else:

    FC_EXONS_PARAMS= " -R -g gene_id -J"




### Target rule input

# Get mutation types to track
MutTypes = config['mut_tracks']
Mutation_Types = MutTypes.split(',')

def get_other_output():

    target = []

    # cB file always gets made
    target.append("results/cB/cB.csv.gz")

    # Tracks always get made
    target.append(expand("results/tracks/{sample}.{mut}.{id}.{strand}.tdf", sample = SAMP_NAMES, mut=Mutation_Types, id=[0,1,2,3,4,5], strand = ['pos', 'min']))

    if config["strategies"]["RSEMp"]:

        target.append("results/transcript_fn/RSEM_plus.csv")

    if config["aligner"] == "star":

        target.append(expand("results/rsem/{SID}.isoforms.results", SID = SAMP_NAMES))

    if config["mutpos"]:

        target.append("results/cB/mutpos.csv.gz")

        target.append("results/cB/mutpos_filtered.csv.gz")

    return target


### Strandedness parameter to pass to mutation counting
if config["strandedness"] == "reverse":

    STRAND = 'R'

else:

    STRAND = 'F'


### Transcript cB helpers

