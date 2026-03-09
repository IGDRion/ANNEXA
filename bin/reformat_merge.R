#! /usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(stringr)

combined <- rtracklayer::import("gffcmp.combined.gtf")
tracking <- read.csv2("gffcmp.tracking", sep = "\t", header = FALSE)

tracking <- tracking %>%
  mutate(
    # Create a combined 'info' column that takes from q1 (V5) if present, else q2 (V6)
    info = coalesce(
      str_extract(V5, "q1:[^\\s]+"),
      str_extract(V6, "q2:[^\\s]+")
    ),
    # Remove q1:/q2: prefix
    info = str_remove(info, "^(q1:|q2:)"),
    # Split the info field by '|' into gene_id and transcript_id
    gene_id = str_extract(info, "^[^|]+"),
    transcript_id = str_extract(info, "(?<=\\|)[^|]+")
  )

meta_combined <- as.data.frame(mcols(combined))
meta_updated <- meta_combined %>%
  left_join(
    tracking %>% 
      select(V1, V2, gene_id_new = gene_id, transcript_id_new = transcript_id),
    by = c("transcript_id" = "V1", "gene_id" = "V2")
  ) %>%
  mutate(
    transcript_id = coalesce(transcript_id_new, transcript_id),
    gene_id = coalesce(gene_id_new, gene_id)
  ) %>%
  select(-transcript_id_new, -gene_id_new)

mcols(combined) <- meta_updated
mcols(combined) <- mcols(combined)[, c("source", "type", "phase", "transcript_id", "gene_id", "exon_number","gene_name")]

rtracklayer::export(combined, "reformated.gtf")