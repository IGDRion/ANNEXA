#! /usr/bin/env Rscript
library(rtracklayer)

#############################################################################
# CLI Args
#############################################################################
args = commandArgs(trailingOnly=TRUE)
class_code_gtf = args[1]
extended_annotation = args[2]

#############################################################################
# Read both GTF and assign gffcompare class_code to extended annotation
#############################################################################
file1 <- readGFF(class_code_gtf)
file2 <- readGFF(extended_annotation)

file2$class_code <- file1$class_code[match(paste(file2$transcript_id, file2$type), paste(file1$transcript_id, file1$type))]

export(file2, extended_annotation)