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
bambu_rec_ndr   <- strsplit(grep('--bambu_rec_ndr*', args, value = TRUE), split = '=')[[1]][[2]]
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
readlist       <- args[8:length(args)]

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
  ),
  type = "message"
)

cat(output, sep = "\n")

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

if (bambu_rec_ndr=="true"){
  # Extract recommended NDR
  ndr_line <- output[grepl("we recommend an NDR of", output)]
  match <- regexec("of ([0-9\\.]+)\\.", ndr_line)
  ndr_value <- as.numeric(regmatches(ndr_line, match)[[1]][2])
  write(ndr_value, file = "rec_ndr.txt")

  # Write filtered GTF and counts
  message("Creating a filtered GTF using a recommended threshold of: ",ndr_value)
  se.filter <- se[mcols(se)$NDR < ndr_value | is.na(mcols(se)$NDR),]
  writeToGTF(rowRanges(se.filter),"extended_annotations.NDR_filter.gtf")
}