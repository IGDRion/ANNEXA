#! /usr/bin/env python3
import argparse
from GTF import GTF

#######################################################
parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
parser.add_argument(
    "--mRNA", help="Path to your mRNA only GTF file.", type=argparse.FileType("r"), required=True
)
parser.add_argument(
    "--lncRNA",
    help="Path to your lncRNA only GTF file.",
    type=argparse.FileType("r"),
    required=True,
)
args = parser.parse_args()

#######################################################
# Parse gtf classed as protein_coding by FEELnc
mRNA_genes = set()

# Add transcript_biotype = protein_coding
for record in GTF.parse_by_line(args.mRNA):
    record.source = "FEELnc"
    record["transcript_biotype"] = "protein_coding"
    record["gene_biotype"] = "protein_coding"

    if record["gene_id"] not in mRNA_genes:
        mRNA_genes.add(record["gene_id"])

    print(record)

#######################################################
# If gene_id also found in protein_coding => protein_coding else lncRNA
for record in GTF.parse_by_line(args.lncRNA):
    record.source = "FEELnc"

    if record["gene_id"] in mRNA_genes:
        record["transcript_biotype"] = "protein_coding"
        record["gene_biotype"] = "protein_coding"
    else:
        record["transcript_biotype"] = "lncRNA"
        record["gene_biotype"] = "lncRNA"

    print(record)
