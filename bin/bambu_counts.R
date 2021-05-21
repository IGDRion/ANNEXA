#! /usr/bin/env Rscript

suppressPackageStartupMessages(library("bambu"))
suppressPackageStartupMessages(library("BSgenome"))

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
args = commandArgs(trailingOnly=TRUE)

output_tag     <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
readCount      <- strsplit(grep('--readCount*', args, value = TRUE), split = '=')[[1]][[2]]
sampleNumber   <- strsplit(grep('--sampleNumber*', args, value = TRUE), split = '=')[[1]][[2]]
ncore          <- strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]]
genomeseq      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
genomeSequence <- Rsamtools::FaFile(genomeseq)
Rsamtools::indexFa(genomeseq)
annot_gtf      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist       <- args[7:length(args)]

print("BAMs:")
readlist

################################################
################################################
## RUN BAMBU                                  ##
################################################
################################################
grlist <- prepareAnnotations(annot_gtf)
se     <- bambu(reads = readlist, 
                annotations = grlist, 
                genome = genomeSequence, 
                ncore = ncore, 
                opt.discovery = list(min.readCount = readCount, 
                                     min.sampleNumber = sampleNumber)
)
writeBambuOutput(se, output_tag)
