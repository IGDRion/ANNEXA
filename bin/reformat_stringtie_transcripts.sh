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

# Step 1: Extract transcript_id, gene_id, and ref_gene_id from reference GTF
awk -F'\t' '$3=="transcript" {
    match($9, /transcript_id "([^"]+)"/, t)
    match($9, /gene_id "([^"]+)"/, g)
    match($9, /ref_gene_id "([^"]+)"/, r)
    if (r[1] != "") {
        print t[1] "\t" r[1]
    } else {
        print t[1] "\t" g[1]
    }
}' "$REFERENCE_GTF" | sort -k1,1 > ref_ids.tmp

# Step 2: Prepare input file (excluding header)
sed '1d' "$INPUT_FILE" | sed 's/,/\t/' | sort -k1,1 > input_sorted.tmp

# Step 3: Join the sorted input file with the reference ids
join -a1 -e "gene_id" -o auto input_sorted.tmp ref_ids.tmp > joined.tmp

# Step 4: Format the final output
(
    # Extract the sample ID from the input file header
    SAMPLE_ID=$(head -n 1 "$INPUT_FILE" | cut -d',' -f2)
    
    # Print the new header
    echo -e "transcript_id\tgene_id\t$SAMPLE_ID"

    # Process and format the rest of the content
    awk '{
        if ($3 == "gene_id") $3 = $2
        print $1 "\t" $3 "\t" $2
    }' joined.tmp
) > "$OUTPUT_FILE"

# Clean up temporary files
rm ref_ids.tmp input_sorted.tmp joined.tmp

echo "Processing complete. Output saved to $OUTPUT_FILE"