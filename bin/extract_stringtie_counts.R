#! /usr/bin/env Rscript
library(IsoformSwitchAnalyzeR)
library(rtracklayer)
library(dplyr)

# CLI Args
args = commandArgs(trailingOnly=TRUE)
av_length = args[1]
ctabs = unlist(strsplit(args[2], " "))
ref = args[3]

stringTieQuant <- importIsoformExpression(
  sampleVector = ctabs,
  readLength = av_length,
  addIsofomIdAsColumn = FALSE
)

myDesign <- data.frame(
  sampleID = colnames(stringTieQuant$abundance),
  condition = gsub('_.*', '', colnames(stringTieQuant$abundance))
)

switchAnalyzeRlist <- importRdata(
  isoformCountMatrix   = stringTieQuant$counts,
  isoformRepExpression = stringTieQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = ref
)

geneCountMatrix <- extractGeneExpression(
  switchAnalyzeRlist,
  extractCounts = TRUE
)

stringTieQuant$counts <- tibble::rownames_to_column(stringTieQuant$counts, var = "transcript_id")
stringTieQuant$counts$gene_id <- with(ref, ifelse(!is.na(ref_gene_id), ref_gene_id, gene_id))[match(stringTieQuant$counts$transcript_id, ref$transcript_id)]
stringTieQuant$counts <- stringTieQuant$counts %>% relocate(gene_id, .after = transcript_id)
write.table(stringTieQuant$counts, "counts_transcript.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(geneCountMatrix, "counts_gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)