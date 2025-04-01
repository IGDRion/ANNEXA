suppressPackageStartupMessages({
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--gtf", help="Path to GTF file")
parser$add_argument("--counts", help="Path to counts file")
parser$add_argument("--output", help="Path to output file (including filename)")
args <- parser$parse_args()

# Sample association
sample_association <- data.frame(
  Sample_Id = c("CFA-MM1", "CFA-MM2", "CFA-MM3", "CFA-MM4", "CFA-HS1", "CFA-HS2", "CFA-HS3", "CFA-OS1", "HSA-MM1", "HSA-MM2"),
  Sample_Name = c("Bear", "CML10", "Popsi", "Twiny", "ChrystalPoumon", "ChrystalRate", "Vico", "Chipie", "HMV2", "WM32")
)

gtf <- rtracklayer::import(args$gtf)
counts <- read.csv2(args$counts, sep = "\t", header = T)
gtf_transcripts <- gtf[gtf$type == "transcript"]

# Create output directory if it doesn't exist
output_dir <- dirname(args$output)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Anonymize sample names in counts
anonymize_sample_name <- function(name) {
  base_name <- sub("_R1.*$", "", name)  # Remove _R1 and anything after it
  base_name <- sub("\\.sorted$", "", base_name)  # Remove .sorted if present
  matched_row <- sample_association[sample_association$Sample_Name == base_name, ]
  if (nrow(matched_row) > 0) {
    return(matched_row$Sample_Id[1])
  } else {
    return(name)  # Keep original name if no match found
  }
}

new_colnames <- sapply(colnames(counts), function(col) {
  if (col %in% c("TXNAME", "GENEID", "transcript_id", "gene_id")) {
    return(col)
  } else {
    return(anonymize_sample_name(col))
  }
})

colnames(counts) <- new_colnames

# Determine the transcript ID column name
tx_id_col <- if("TXNAME" %in% colnames(counts)) "TXNAME" else "transcript_id"

# Extract data ------------------------------------------------------------

# Function to filter gtf_transcripts for a given sample
filter_transcripts <- function(sample_name, counts_data, gtf_data, tx_id_col) {
  expressed_transcripts <- counts_data[counts_data[[sample_name]] > 0, tx_id_col]
  filtered_gtf <- gtf_data[gtf_data$transcript_id %in% expressed_transcripts]
  return(filtered_gtf)
}

# Get sample names (excluding transcript and gene ID columns)
sample_names <- setdiff(colnames(counts), c("TXNAME", "GENEID", "transcript_id", "gene_id"))

# Create a list to store filtered gtf objects for each sample
filtered_gtf_list <- list()

# Loop through each sample and create filtered gtf objects
for (sample in sample_names) {
  filtered_gtf_list[[sample]] <- filter_transcripts(sample, counts, gtf_transcripts, tx_id_col)
  print(paste("Created filtered gtf for", sample, "with", length(filtered_gtf_list[[sample]]), "transcripts"))
}

# Plot data ---------------------------------------------------------------

# Function to calculate proportions of class codes
get_class_code_proportions <- function(gtf) {
  gtf %>%
    as.data.frame() %>%
    filter(class_code %in% c("=", "x", "u", "j", "k")) %>%
    count(class_code) %>%
    mutate(proportion = n / sum(n))
}

# Calculate proportions for each sample
plot_data <- lapply(names(filtered_gtf_list), function(sample) {
  props <- get_class_code_proportions(filtered_gtf_list[[sample]])
  props$sample <- sample
  return(props)
}) %>% bind_rows()

# Create the stacked bar plot
plot <- ggplot(plot_data, aes(x = sample, y = proportion, fill = class_code)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Sample", y = "Proportion", fill = "Class Code") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent)

# Save the plot to the specified output file
pdf(args$output)
print(plot)
dev.off()
message(paste("Plot saved successfully to", args$output))
