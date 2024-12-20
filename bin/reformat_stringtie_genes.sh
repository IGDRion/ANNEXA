#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file> <reference_gtf> <output_file>"
    exit 1
fi

# Assign command-line arguments to variables
INPUT_FILE="$1"
REFERENCE_GTF="$2"
OUTPUT_FILE="$3"

# Check if input files exist
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

if [ ! -f "$REFERENCE_GTF" ]; then
    echo "Error: Reference GTF file $REFERENCE_GTF does not exist."
    exit 1
fi

# Step 1: Extract gene_id, gene_name, and ref_gene_id mappings from reference GTF
awk -F'\t' '$3=="transcript" {
    match($9, /gene_id "([^"]+)"/, g)
    match($9, /gene_name "([^"]+)"/, n)
    match($9, /ref_gene_id "([^"]+)"/, r)
    if (g[1] != "" && n[1] != "" && r[1] != "") 
        print g[1] "\t" n[1] "\t" r[1]
}' "$REFERENCE_GTF" | sort -u > gene_mappings.tmp

# Step 2: Process the input file
(
    # Dynamically extract and print the header from the input file
    head -n 1 "$INPUT_FILE" | sed 's/,/\t/'

    # Process the rest of the file
    awk -F',' 'NR>1 {
        split($1, a, "|")
        if (a[2] != "") {
            if (a[1] ~ /^ENSCAFG/) {
                print a[1] "\t" $2
            } else {
                cmd = "awk -v gid=\"" a[1] "\" -v gname=\"" a[2] "\" '\''$1 == gid && $2 == gname {print $3}'\'' gene_mappings.tmp"
                cmd | getline ref_gene_id
                close(cmd)
                if (ref_gene_id != "") {
                    print ref_gene_id "\t" $2
                } else {
                    print $1 "\t" $2  # Keep original if no mapping found
                }
            }
        } else {
            print $1 "\t" $2  # Keep original for MSTRG.* without |
        }
    }' "$INPUT_FILE"
) > "$OUTPUT_FILE"

# Clean up
rm gene_mappings.tmp