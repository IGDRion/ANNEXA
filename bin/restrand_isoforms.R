#! /usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(stringr)

#############################################
# CLI Args
args = commandArgs(trailingOnly=TRUE)
gtf_file = args[1]
tx_tool = args[2]
output = args[3]

# Define prefix to identify novel transcripts
if (tx_tool == "stringtie2") {
  prefix = "MSTRG"
} else{
    prefix = "Bambu"
}

gtf <- rtracklayer::readGFF(gtf_file)

# Restrand isoforms of known genes ----------------------------------------

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

# Change strand of novel isoforms in gtf dataframe
gtf <- gtf %>%
  mutate(strand = case_when(
    transcript_id %in% to_change$transcript_id ~ to_change$strand[match(transcript_id, to_change$transcript_id)],
    TRUE ~ strand
  ))

# Restrand novel tx of novel genes ----------------------------------------

# Look for novel tx of novel genes
novel <- gtf %>% filter(str_detect(transcript_id, prefix) & str_detect(gene_id, prefix))

# Restrand novel transcripts
# If all isoforms are unstraded, keep unstranded
# If some are stranded, use strand as new strand
# If all three strands are there, use majority rule
new_gene_strands <- novel %>%
  group_by(gene_id) %>%
  summarize(gene_strand = if(all(strand == "*")) "*" 
            else {
              strands <- table(strand[strand != "*"])
              if(length(strands) > 1 && max(strands) == min(strands)) "*"
              else names(which.max(strands))
            })

# Change strand in gtf
gtf <- gtf %>%
  left_join(new_gene_strands, by="gene_id") %>%
  mutate(strand = coalesce(gene_strand, strand)) %>%
  select(-gene_strand)

# Export gtf to new restranded file
export(gtf, output)