#! /usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(stringr)

# This script reassigns the strand of novel isoforms identified by bambu or stringtie.
# For novel isoforms of known reference genes, it assigns the strand of the reference gene \
# to its novel isoform.
# Example:
# chr21	StringTie	transcript	8380666	8451371	1000	+	.	gene_id "ENSG00000280441.3"; transcript_id "ENST00000652379.1"; gene_name "ENSG00000280441"; ref_gene_id "ENSG00000280441.3"; 
# chr21	StringTie	exon	8380666	8380972	1000	+	.	gene_id "ENSG00000280441.3"; transcript_id "ENST00000652379.1"; exon_number "1"; gene_name "ENSG00000280441"; ref_gene_id "ENSG00000280441.3"; 
# chr21	StringTie	transcript	8396906	8397356	1000	.	.	gene_id "ENSG00000280441.3"; transcript_id "MSTRG.19526.1"; 
# chr21	StringTie	exon	8396906	8397356	1000	.	.	gene_id "ENSG00000280441.3"; transcript_id "MSTRG.19526.1"; exon_number "1"; 
# chr21	StringTie	transcript	8403347	8404298	1000	-	.	gene_id "ENSG00000280441.3"; transcript_id "MSTRG.19528.1"; 
# chr21	StringTie	exon	8403347	8404298	1000	-	.	gene_id "ENSG00000280441.3"; transcript_id "MSTRG.19528.1"; exon_number "1"; 
# The two novel transcripts (identified by MSTRG. or BambuTx) will be restranded to the + strand based on the strand of the ref_gene_id.

#############################################
# CLI Args
args = commandArgs(trailingOnly=TRUE)
gtf_file = args[1]
tx_tool = args[2]
output = args[3]

# Define prefix to identify novel transcripts
if (tx_tool == "stringtie") {
  prefix = "MSTRG"
} else{
    prefix = "Bambu"
}

gtf <- rtracklayer::readGFF(gtf_file)

# Check if gtf already has ref_gene_id feature (Stringtie), if not (Bambu), create it based on gene_id column
if (!"ref_gene_id" %in% colnames(gtf)){
  gtf$ref_gene_id <- ifelse(!grepl("Bambu|MSTRG", gtf$gene_id), 
                                  gtf$gene_id, 
                                  NA)
}

# Get strand of reference genes
ref_strand <- unique(gtf %>% filter(!is.na(ref_gene_id)) %>% select(ref_gene_id, strand))
# Find novel isoforms
to_change <- gtf %>% filter(str_detect(transcript_id, prefix) & !str_detect(gene_id, prefix))
# Match and apply new strand (if different) to novel isoforms
matches <- match(to_change$gene_id, ref_strand$ref_gene_id)
to_change$strand <- ref_strand$strand[matches]

# Change strand of novel isoforms in gtf dataframe and export as new gtf
gtf <- gtf %>%
  mutate(strand = case_when(
    transcript_id %in% to_change$transcript_id ~ to_change$strand[match(transcript_id, to_change$transcript_id)],
    TRUE ~ strand
  ))

export(gtf, output)