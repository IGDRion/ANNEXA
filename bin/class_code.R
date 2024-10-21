#! /usr/bin/env Rscript
library(rtracklayer)
library(dplyr)

#############################################################################
# CLI Args
#############################################################################
args = commandArgs(trailingOnly=TRUE)
class_code_gtf = args[1]
extended_annotation = args[2]
output = args[3]

#############################################################################
# Read both GTF and assign gffcompare class_code to extended annotation
#############################################################################
cc_gtf <- readGFF(class_code_gtf)
ext_anno <- readGFF(extended_annotation)

# Add class_code
ext_anno$class_code <- cc_gtf$class_code[match(paste(ext_anno$transcript_id, ext_anno$type), paste(cc_gtf$transcript_id, cc_gtf$type))]

# Change order of transcript_id column
ext_anno <- ext_anno %>% relocate(transcript_id, .after = gene_id)

export(ext_anno, output)