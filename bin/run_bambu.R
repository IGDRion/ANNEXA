#! /usr/bin/env Rscript

suppressPackageStartupMessages(library("bambu"))
suppressPackageStartupMessages(library("BSgenome"))

################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
args = commandArgs(trailingOnly = TRUE)

output_tag     <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore          <- strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]]
genomeseq      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
genomeSequence <- Rsamtools::FaFile(genomeseq)
Rsamtools::indexFa(genomeseq)
annot_gtf      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist       <- args[5:length(args)]

print("BAMs:")
readlist

################################################
## RUN BAMBU                                  ##
################################################
grlist <- prepareAnnotations(annot_gtf)
se     <- bambu(
  reads = readlist,
  annotations = grlist,
  genome = genomeSequence,
  ncore = ncore,
  verbose = TRUE,
  NDR = 1,
  opt.discovery = list(min.txScore.singleExon = 0)
)

# Extract NDR
tx_infos   <- se@rowRanges@elementMetadata
new_tx_idx <- tx_infos[tx_infos$novelTranscript == "TRUE" & tx_infos$txClassDescription != "annotation", c(1,2,3) ]
write.csv(
  new_tx_idx,
  "bambu_ndr.csv",
  quote = FALSE,
  row.names = FALSE
)

# Write GTF and counts
writeBambuOutput(se, output_tag)
