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
bambu_strand   <- strsplit(grep('--bambu_strand*', args, value = TRUE), split = '=')[[1]][[2]]
bambu_singleexon   <- strsplit(grep('--bambu_singleexon*', args, value = TRUE), split = '=')[[1]][[2]]
if (bambu_strand=="true"){
  bambu_strand=TRUE
} else {
  bambu_strand=FALSE
}
if (bambu_singleexon=="true"){
  bambu_singleexon=TRUE
} else {
  bambu_singleexon=FALSE
}
readlist       <- args[7:length(args)]

print("BAMs:")
readlist

################################################
## RUN BAMBU                                  ##
################################################
grlist <- prepareAnnotations(annot_gtf)

if (bambu_singleexon==TRUE) {
  opt_discovery <- list(min.txScore.singleExon = 0)
} else {
  opt_discovery <- NULL
}

output <- capture.output(
  se     <- bambu(
    reads = readlist,
    annotations = grlist,
    genome = genomeSequence,
    ncore = ncore,
    verbose = TRUE,
    NDR = 1,
    stranded = bambu_strand,
    opt.discovery = opt_discovery
  )
)

cat(output, sep = "\n")

# Extract recommended NDR
ndr_line <- output[grepl("we recommend an NDR of", output)]
ndr_value <- as.numeric(gsub(".*we recommend an NDR of ([0-9.]+).*", "\\1", ndr_line))
write(ndr_value, file = "rec_ndr.txt")

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
