
### 





### PROCESSING PARAMETERS PASSED BY USER

str = "--outSAMattributes NH HI AS jI --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts"

# Add necessary tags and change other settings as necessary
args = str.split()

# Tags to always include
    # Technically, NH, HI, AS, and nM are added by default
    # NH = number of places the read aligns
    # HI = Numerical index that is useful to separate multi-mapping reads
    # NM = Edit distance between read and reference
    # AS = Alignment score
    # MD = Used for mutation counting
    # nM = Number of mismatches
sam_attributes = set(["NH", "HI", "AS", "NM", "MD"])

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

# Process --quantMode TranscriptomeSAM
quant_mode_str = "--quantMode TranscriptomeSAM"
if quant_mode_str not in ' '.join(args):
    args.append("--quantMode")
    args.append("TranscriptomeSAM")
    

# Example usage
input_string = "--outSAMattributes NH AS --otherFlag value"
modified_string = modify_star_arguments(input_string)
print(modified_string)




###  my crack at it

STAR_PARAMS = "--outSAMattributes NH HI AS jI --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"

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
sam_attributes = set(["NH", "HI", "AS", "NM", "MD"])

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
if "--outSAMType" in args:
    
    index = args.index("--outSAMType")

    # Replace the existing outSAMType values with the desired ones
    args[index + 1:index + 3] = ["BAM", "SortedByCoordinate"]

else:
    
    args.extend(["--outSAMType", "BAM", "SortedByCoordinate"])


' '.join(args)