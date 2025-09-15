library(dplyr)
library(ggplot2)
library(ggtranscript)
library(rtracklayer)
library(tidyr)

# Import data -------------------------------------------------------------

# Import gtf
ref_gtf <- rtracklayer::import("Projets/ggtranscript/annotation_19.gtf")
new_gtf <- rtracklayer::import("Projets/ggtranscript/extended_annotations.filter.gtf")

# Find overlaps for novel exons
new_gtf$diff_type <- ifelse(new_gtf %in% ref_gtf, "in_ref", "not_in_ref")

# Convert to tibble
ref_gtf <- ref_gtf %>% dplyr::as_tibble()
ref <- ref_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_id,
    transcript_id,
    transcript_biotype
  )

new_gtf <- new_gtf %>% dplyr::as_tibble()
new <- new_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_id,
    transcript_id,
    transcript_biotype,
    diff_type
  )

# Create diffs
exons <- new %>% filter(type == "exon")
diffs <- as.data.frame(list(
  seqnames = exons$seqnames,
  start = exons$start,
  end = exons$end,
  width = abs(exons$end - exons$start + 1),
  strand = exons$strand,
  type = "diff",
  diff_type = exons$diff_type,
  transcript_id = exons$transcript_id
))

# Import junction data and match ------------------------------------------

#Import junction data
junction_data <- read.csv2("Projets/ggtranscript/test.jcounts", header = T, sep = "\t")
junctions <- data.frame(
  seqnames = junction_data$Site1_chr,
  start    = junction_data$Site1_location,
  end      = junction_data$Site2_location,
  strand   = junction_data$Site1_strand,
  mean_count = rowMeans(junction_data[, grep("bam", colnames(junction_data))])
)

exons <- ref %>% dplyr::filter(type == "exon")
sorted <- exons %>%
  arrange(transcript_id, start) %>%
  group_by(transcript_id)

exon_boundaries <- sorted %>%
  mutate(next_start = lead(start), next_transcript = lead(transcript_id)) %>%
  filter(transcript_id == next_transcript) %>%
  select(seqnames, exon_end = end, exon_next_start = next_start, transcript_id)

# Match gene_id
### Match gene_id 
transcript_to_gene <- ref %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id) %>%
  distinct()

find_transcripts_for_junction <- function(junc_start, junc_end, junc_seqname, boundaries_df) {
  transcripts <- boundaries_df %>%
    filter(seqnames == junc_seqname,
           exon_end == junc_start,
           exon_next_start == junc_end) %>%
    pull(transcript_id) %>%
    unique()
  if(length(transcripts) == 0) {
    return(NA_character_)
  } else {
    return(paste(transcripts, collapse = ","))
  }
}

# Annotate junctions with transcript names, duplicate if multiple isoforms for a gene, add gene_id
annotated_junctions <- junctions %>%
  rowwise() %>%
  mutate(transcript_id = find_transcripts_for_junction(start, end, seqnames, exon_boundaries)) %>%
  ungroup() %>%
  separate_rows(transcript_id, sep = ",") %>%
  mutate(transcript_id = trimws(transcript_id)) %>%
  left_join(transcript_to_gene, by = "transcript_id")

# Plot --------------------------------------------------------------------
wanted = "ENSCAFG00000004151"
example <- new %>% dplyr::filter(gene_id == wanted)
example_exons <-  example %>% filter(type =="exon")
example_cds <- example %>% filter(type == "CDS")

example_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    data = example_cds,
    aes(fill = diff_type),
    linewidth = 0.5
  ) +
  geom_intron(
    data = to_intron(example_exons, "transcript_id"),
    aes(strand = strand)
  ) + 
  geom_junction(
    data = annotated_junctions %>% filter(gene_id == wanted),
    aes(size = mean_count),
    junction.y.max = 0.5
  ) +
  geom_junction_label_repel(
    data = annotated_junctions %>% filter(gene_id == wanted),
    aes(label = round(mean_count, 2)),
    junction.y.max = 0.5
  ) + 
  scale_size_continuous(range = c(0.1, 1)) +
  ggtitle(paste("Gene id:"), wanted) +
  geom_range(
    data = diffs %>% filter(transcript_id %in% example$transcript_id),
    aes(fill = diff_type),
    alpha = 0.2
  ) +
  geom_text(
    data = add_exon_number(example_exons, "transcript_id"),
    aes(
      x = (start + end) / 2, # plot label at midpoint of exon
      label = exon_number
    ),
    size = 6,
    nudge_y = 0.4
  ) + 
  theme_bw() + 
  scale_x_continuous(name = "Position") + 
  scale_y_discrete(name = "Transcript name")
