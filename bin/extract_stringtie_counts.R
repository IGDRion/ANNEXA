#! /usr/bin/env Rscript
suppressMessages(suppressWarnings(library(IsoformSwitchAnalyzeR)))
suppressMessages(suppressWarnings(library(rtracklayer)))
library(dplyr)
library(tibble)

# CLI Args
args = commandArgs(trailingOnly=TRUE)

av_length <- as.numeric(args[which(args == "--av_length") + 1])
ref <- args[which(args == "--ref") + 1]
ctabs <- args[(which(args == "--ctabs") + 1):(which(args == "--ref") - 1)]

# -------------------------------------------------------------------------
# Extracting counts with IsoformSwitchAnalyzeR
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

# -------------------------------------------------------------------------
# Renaming Stringtie MSTRG gene_id to original reference gene_id, if possible, and
# add unexpressed genes (left out by IsoformSwitchAnalyzeR) back to geneCountMatrix

ref <- rtracklayer::import(ref)

## Rename missing MSTRG to ref_gene_id in geneCountMatrix when unambiguous
## Ambiguous gene_ids are left as they are (with MSTRG id) and referenced in a new list
# Rename missing MSTRG to ref_gene_id in geneCountMatrix
mstrg_mapping <- as.data.frame(ref[startsWith(ref$gene_id, "MSTRG") & !is.na(ref$ref_gene_id), c("gene_id", "ref_gene_id")])
ambiguous_mappings <- mstrg_mapping %>%
  group_by(gene_id) %>%
  summarize(ref_gene_ids = list(unique(ref_gene_id)), count = n_distinct(ref_gene_id)) %>%
  filter(count > 1)
unambiguous_mappings <- mstrg_mapping %>%
  anti_join(ambiguous_mappings, by = "gene_id")
mstrg_lookup <- setNames(unambiguous_mappings$ref_gene_id, unambiguous_mappings$gene_id)
geneCountMatrix$gene_id <- ifelse(geneCountMatrix$gene_id %in% names(mstrg_lookup),
                                  mstrg_lookup[geneCountMatrix$gene_id],
                                  geneCountMatrix$gene_id)

ambiguous_list <- setNames(
  lapply(ambiguous_mappings$ref_gene_ids, function(x) paste(x, collapse = ", ")),
  ambiguous_mappings$gene_id
)
ambiguous_df <- data.frame(
  gene_id = names(ambiguous_list),
  ref_gene_id = unlist(ambiguous_list),
  stringsAsFactors = FALSE
)
if (nrow(ambiguous_mappings) > 0){
write.table(ambiguous_df, file = "ambiguous_mappings.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  warning("Some gene_ids could not be mapped to a ref_gene_id because they are assigned by Stringtie to multiple reference genes.\
  These genes are listed in 'ambiguous_mappings.tsv'")
}

# Search missing_genes, and genes already fixed
missing_genes <- setdiff(unique(ref$gene_id), unique(geneCountMatrix$gene_id))
fixed_genes <- setdiff(unique(geneCountMatrix$gene_id), unique(ref$gene_id))

# Convert fixed_genes (ref_id) back to gene_id names and remove them from missing_genes
gene_id_df <- data.frame(gene_id = ref$gene_id, ref_gene_id = ref$ref_gene_id)
gene_id_df <- unique(gene_id_df[!is.na(gene_id_df$ref_gene_id), ])
gene_id_map <- setNames(gene_id_df$gene_id, gene_id_df$ref_gene_id)
converted_genes <- gene_id_map[fixed_genes]
converted_genes <- unname(converted_genes[!is.na(converted_genes)])
missing_genes <- setdiff(missing_genes, converted_genes)

# Convert missing_genes to ref_gene_id, if there is one
gene_map <- data.frame(gene_id = ref$gene_id, ref_gene_id = ref$ref_gene_id)
gene_map <- unique(gene_map[!is.na(gene_map$gene_id) & !is.na(gene_map$ref_gene_id), ])
id_map <- setNames(gene_map$ref_gene_id, gene_map$gene_id)
renamed_missing_genes <- ifelse(
  missing_genes %in% names(id_map),  # Check if missing_genes exist in id_map
  id_map[missing_genes],            # Replace with corresponding ref_gene_id
  missing_genes                     # Keep original value if no mapping exists
)
renamed_missing_genes <- unname(renamed_missing_genes)

# Add missing genes to geneCountMatrix
missing_df <- data.frame(
  gene_id = renamed_missing_genes,
  gene_name = "NA",
  matrix(0, nrow = length(renamed_missing_genes), ncol = ncol(geneCountMatrix) - 2)
)
colnames(missing_df)[3:ncol(missing_df)] <- colnames(geneCountMatrix)[3:ncol(geneCountMatrix)]
updated_geneCountMatrix <- rbind(geneCountMatrix, missing_df)
updated_geneCountMatrix <- updated_geneCountMatrix[order(updated_geneCountMatrix$gene_id), ]
rownames(updated_geneCountMatrix) <- NULL
updated_geneCountMatrix <- updated_geneCountMatrix %>% select(gene_id, gene_name, sort(setdiff(names(.), c("gene_id", "gene_name"))))

# -------------------------------------------------------------------------
# Add gene_id to stringtieQuant$counts (transcript count)
lookup <- as.data.frame(ref) %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id, ref_gene_id) %>%
  distinct() %>%
  group_by(gene_id) %>%
  mutate(
    unique_ref_gene_id = n_distinct(ref_gene_id[!is.na(ref_gene_id)]) == 1,
    common_ref_gene_id = if(unique_ref_gene_id[1] && any(!is.na(ref_gene_id))) 
      ref_gene_id[!is.na(ref_gene_id)][1] 
    else 
      NA_character_
  ) %>%
  ungroup()

updated_stringTieQuant <- stringTieQuant$counts %>%
  rownames_to_column("transcript_id") %>%
  left_join(lookup, by = "transcript_id") %>%
  mutate(
    final_gene_id = case_when(
      !is.na(common_ref_gene_id) ~ common_ref_gene_id,
      !is.na(ref_gene_id) ~ ref_gene_id,
      TRUE ~ gene_id
    )
  ) %>%
  select(transcript_id, names(stringTieQuant$counts), final_gene_id) %>%
  rename(gene_id = final_gene_id) %>%
  column_to_rownames("transcript_id")
updated_stringTieQuant <- tibble::rownames_to_column(updated_stringTieQuant, var = "transcript_id")
updated_stringTieQuant <- updated_stringTieQuant %>% relocate(gene_id, .after = transcript_id)
updated_stringTieQuant <- updated_stringTieQuant %>% select(transcript_id, gene_id, sort(setdiff(names(.), c("transcript_id", "gene_id"))))

# -------------------------------------------------------------------------
write.table(updated_stringTieQuant, "counts_transcript.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(updated_geneCountMatrix[,-2], "counts_gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)

message("Finished renaming process. Gene and transcripts counts are written in 'counts_gene.txt' and 'counts_transcript.txt'")