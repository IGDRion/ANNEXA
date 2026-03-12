#! /usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(stringr)

combined <- rtracklayer::import("gffcmp.combined.gtf")
tracking <- read.csv2("gffcmp.tracking", sep = "\t", header = FALSE)

tracking <- tracking %>%
  mutate(
    q1 = str_extract(V5, "q1:[^\\s]+"),
    q2 = str_extract(V6, "q2:[^\\s]+"),
    
    q1 = str_remove(q1, "^q1:"),
    q2 = str_remove(q2, "^q2:"),
    
    gene_id_q1 = ifelse(!is.na(q1), str_extract(q1, "^[^|]+"), NA),
    transcript_id_q1 = ifelse(!is.na(q1), str_extract(q1, "(?<=\\|)[^|]+"), NA),
    gene_id_q2 = ifelse(!is.na(q2), str_extract(q2, "^[^|]+"), NA),
    transcript_id_q2 = ifelse(!is.na(q2), str_extract(q2, "(?<=\\|)[^|]+"), NA),
    
    gene_id = coalesce(gene_id_q1, gene_id_q2),
    transcript_id = coalesce(transcript_id_q1, transcript_id_q2),
    
    alt = ifelse(
      !is.na(q2) & q2 != "-" &
        grepl("^(BambuTx|MSTRG\\.)", transcript_id_q2, perl = TRUE),
      transcript_id_q2,
      NA_character_
    )
  )

meta_combined <- as.data.frame(mcols(combined))
meta_updated <- meta_combined %>%
  left_join(
    tracking %>% 
      select(V1, V2, gene_id_new = gene_id, transcript_id_new = transcript_id, alt),
    by = c("transcript_id" = "V1", "gene_id" = "V2")
  ) %>%
  mutate(
    transcript_id = coalesce(transcript_id_new, transcript_id),
    gene_id = coalesce(gene_id_new, gene_id)
  ) %>%
  select(-transcript_id_new, -gene_id_new)

mcols(combined) <- meta_updated
mcols(combined) <- mcols(combined)[, c("source", "type", "phase", "transcript_id", "gene_id", "exon_number","gene_name","alt")]

rtracklayer::export(combined, "reformated.gtf")