#! /usr/bin/env Rscript
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggridges)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(viridis)
library(ggpattern)
library(stringr)

brew = c("#008387", "#a0d3db", "#ad8500", "#d9cb98")
palette = c("#00AFBB", "#FC4E07", "#E7B800")

theme_set(
  theme_pubr(base_size = 15) +
    theme(
      panel.grid.major = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      strip.text = element_text(face = "bold"),
      plot.margin = unit(c(1, 1, 1, 0.5), "cm"),
      plot.title = element_text(
        color = "black",
        size = 18,
        face = "bold",
        hjust = 0.5
      )
    )
)

#############################################################################
# CLI Args
#############################################################################
args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
ANNEXA_version = args[2]
command = args[3]

#############################################################################
# GENE level
#############################################################################
gene <- read.csv(paste0(prefix,".gene.stats"), header = TRUE)
gene$gene_biotype[gene$gene_biotype == "protein_coding"] <- "mRNA"

# Length
med_length = gene %>%
  group_by(discovery, gene_biotype) %>%
  summarize(median = median(length), gene_biotype = gene_biotype)

len = gene %>%
  ggplot(aes(x = length, fill = paste(gene_biotype, discovery))) +
  ggtitle("Gene length distribution") +
  geom_vline(data = med_length,
             aes(
               xintercept = median,
               color = paste(gene_biotype, discovery)
             ),
             size = 1) +
  geom_density(alpha = 0.8) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  facet_grid(gene_biotype ~ .) +
  scale_fill_manual(values = brew) +
  scale_color_manual(values = brew) +
  guides(color = FALSE) +
  guides(fill = guide_legend("Source")) +
  xlab("Gene length (nt)")

# Nb isoformes
iso = gene %>%
  mutate(isoformes = ifelse(nb_transcripts == 1, "1", "2+")) %>%
  ggplot(aes(x = discovery, fill = interaction(discovery,gene_biotype), pattern = isoformes)) +
  ggtitle("Proportion of mono versus multi-isoform genes") +
  geom_bar_pattern(position = "fill",
           alpha = 0.8,
           colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~ gene_biotype) +
  scale_fill_manual(values = brew) +
  scale_pattern_manual(values = c("1" = "stripe", "2+" = "none")) +
  guides(fill = "none", pattern = guide_legend(title = "Isoforms")) +
  theme(axis.title.x = element_blank())

# Nb  of genes per biotype and statuts known/novel
count = ggplot(data = gene, aes(x = gene_biotype, fill = paste(gene_biotype, discovery))) +
  ggtitle("Number of genes") +
  geom_bar(colour = "black",
           alpha = 0.8,
           position = "dodge") +
  scale_fill_manual(values = brew) +
  guides(fill = guide_legend("Source")) +
  geom_text(
    stat = 'count',
    aes(label = ..count..),
    vjust = -0.2,
    position = position_dodge(width = 0.9)
  ) +
  ylab("Number of gene") +
  theme(axis.title.x = element_blank())

# Distribution of gene counts in ridge plot
gene_counts = gene %>%
  filter(validate_by >= 5) %>%
  ggplot() +
  ggtitle("Distribution of gene counts") +
  geom_density_ridges2(aes(
    x = validate_by,
    y = discovery,
    fill = paste(gene_biotype, discovery)
  ), alpha = 0.7) +
  theme(legend.position = "top") +
  facet_grid(gene_biotype ~ .) +
  scale_fill_manual(values = brew) +
  guides(fill = guide_legend("Source")) +
  scale_x_log10() +
  xlab("Gene counts") +
  theme(axis.title.y = element_blank()) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  )

# Hist of gene validate by X samples
gene_validate = gene %>%
  filter(presents_in_sample > 0) %>%
  ggplot(aes(presents_in_sample, fill = paste(gene_biotype, discovery))) +
  ggtitle("Number of genes according to his\nbreadth of expression (gene count >1)") +
  geom_bar(color = "black", position = position_dodge()) +
  facet_grid(gene_biotype ~ .) +
  scale_y_log10() +
  scale_fill_manual(values = brew) +
  guides(fill = guide_legend("Source")) +
  xlab("Number of samples") +
  ylab("Number of genes") +
  geom_text(
    stat = 'count',
    aes(label = ..count..),
    vjust = -0.2,
    position = position_dodge(width = 0.9)
  )

# Distribution counts by sample
gene_counts_sample = gene %>%
  filter(presents_in_sample > 0 & validate_by >= 5) %>%
  ggplot(aes(
    x = validate_by,
    y = factor(presents_in_sample),
    fill = paste(gene_biotype, discovery)
  )) +
  ggtitle("Gene counts separated by number of samples\nwith at least 1 count") +
  geom_density_ridges2(alpha = 0.6) +
  scale_x_log10() +
  facet_grid(gene_biotype ~ .) +
  scale_fill_manual(values = brew) +
  scale_color_manual(values = brew) +
  guides(color = FALSE) +
  guides(fill = guide_legend("Source")) +
  xlab("Gene counts") +
  ylab("Number of samples") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  )

# 5'-3' extensions count
gene_ext = gene %>%
  filter(ext_5 > 0 | ext_3 > 0) %>%
  mutate(ext = ifelse(ext_5 > 0 & ext_3 > 0, "5'-3'",
                      ifelse(ext_5 > 0, "5'", "3'"))) %>%
  ggplot(aes(x = gene_biotype, fill = gene_biotype, pattern = ext, pattern_angle = ext)) +
  ggtitle("Number of 5' and 3' gene extensions") +
  geom_bar_pattern(colour = "black",
           alpha = 0.8,
           position = position_dodge()) +
  scale_fill_manual(values = brew[c(2,4)]) +
  scale_pattern_manual(values = c("stripe","stripe","crosshatch")) +
  scale_pattern_angle_manual(values = c(45, -45, 45)) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_blank()) +
  labs(y = "Count") +
  guides(fill = "none",
         pattern = guide_legend(title = "Gene extension"),
         pattern_angle = guide_legend(title = "Gene extension")) +
  scale_y_log10() +
  annotation_logticks(sides = "l")

# 5'-3' extensions distribution
labels = c("5'", "3'")
names(labels) = c("ext_5", "ext_3")
gene_ext_dist = gene %>%
  filter(ext_5 > 0 | ext_3 > 0) %>%
  select(c("ext_5", "ext_3", "gene_biotype", "discovery")) %>%
  melt() %>%
  ggplot(aes(x = value, fill = paste(gene_biotype, discovery))) +
  ggtitle("Distribution of 5' and 3' gene extensions\n(at genomic level)") +
  geom_density(alpha = 0.5) +
  facet_grid(variable ~ ., labeller = labeller(variable = labels)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  theme(legend.position = "top") +
  guides(fill = guide_legend("Source")) +
  xlab("Length of gene extension (nt)") + ylab("Density") +
  scale_fill_manual(values = c(brew[1], brew[3]))


#############################################################################
# TRANSCRIPT
#############################################################################
transcript <- read.csv(paste0(prefix, ".transcript.stats"), header = TRUE)
lncRNA_biotypes <- c("lncRNA",
                    "antisense",
                    "non-coding",
                    "lnc_RNA",
                    "ncRNA")
mRNA_biotypes <- c("protein_coding", "mRNA")

transcript <- transcript %>%
  mutate(transcript_biotype = if_else(
    transcript_biotype %in% lncRNA_biotypes,
    "lncRNA",
    if_else(
      transcript_biotype %in% mRNA_biotypes,
      "mRNA",
      transcript_biotype
    )
  )) %>%
  filter(transcript_biotype %in% c("mRNA", "lncRNA"))

dim(transcript)
# Length distrib
med_length_tx <- transcript %>%
  group_by(tx_discovery, transcript_biotype) %>%
  summarize(median = median(length), transcript_biotype = transcript_biotype)

tx_len = transcript %>%
  ggplot(aes(
    x = length,
    fill = paste(transcript_biotype, tx_discovery)
  )) +
  ggtitle("Transcript length distribution") +
  geom_vline(data = med_length_tx,
             aes(
               xintercept = median,
               color = paste(transcript_biotype, tx_discovery)
             ),
             size = 1) +
  geom_density(alpha = 0.7) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  facet_grid(transcript_biotype ~ .) +
  guides(fill = guide_legend("Source")) +
  scale_fill_manual(values = brew) +
  scale_color_manual(values = brew) +
  guides(color = FALSE) +
  xlab("Transcript length (nt)")

# Mono vs multiexonique
tx_ex = transcript %>%
  mutate(exons = ifelse(nb_exons == 1, "1", "2+")) %>%
  ggplot(aes(x = tx_discovery, fill = interaction(tx_discovery, transcript_biotype), pattern = exons)) +
  ggtitle("Proportion of mono versus multi exonic transcripts") +
  geom_bar_pattern(position = "fill",
           alpha = 0.8,
           colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~ transcript_biotype) +
  scale_fill_manual(values = brew) +
  scale_pattern_manual(values = c("1" = "stripe", "2+" = "none")) +
  guides(fill = "none", pattern = guide_legend(title = "Exons")) +
  theme(axis.title.x = element_blank())

# Count
tx_count = transcript %>%
  ggplot(aes(
    x = transcript_biotype,
    fill = paste(transcript_biotype, tx_discovery, "in", gene_discovery, "gene")
  )) +
  ggtitle("Number of transcripts") +
  geom_bar(colour = "black",
           alpha = 0.8,
           position = "dodge") +
  scale_fill_manual(
    values = c("#006063", brew[1], brew[2], "#876700", brew[3], brew[4]),
    labels = c(
      "K in K gene",
      "N in K gene",
      "N in N gene",
      "K in K gene",
      "N in K gene",
      "N in N gene"
    )
  ) +
  guides(fill = guide_legend("K=known, N=novel")) +
  geom_text(
    stat = 'count',
    aes(label = ..count..),
    vjust = -0.2,
    position = position_dodge(width = 0.9)
  ) +
  theme(legend.title = element_text(size = 10), axis.title.x = element_blank()) +
  ylab("Number of transcripts")

#############################################################################
# EXON
#############################################################################
exon = read.csv(paste0(prefix,".exon.stats"), header = T)
exon = exon %>%
  mutate(exon_biotype = if_else(
    exon_biotype %in% lncRNA_biotypes,
    "lncRNA",
    if_else(
      exon_biotype %in% mRNA_biotypes,
      "mRNA",
      exon_biotype
    )
  )) %>%
  filter(exon_biotype %in% c("mRNA", "lncRNA"))

# Length distrib
med_length_ex = exon %>%
  group_by(discovery, exon_biotype) %>%
  summarize(median = median(length), exon_biotype = exon_biotype)

ex_len = exon %>%
  ggplot(aes(x = length, fill = paste(exon_biotype, discovery))) +
  ggtitle("Exon length distribution") +
  geom_vline(data = med_length_ex,
             aes(
               xintercept = median,
               color = paste(exon_biotype, discovery)
             ),
             size = 1) +
  geom_density(alpha = 0.7) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  facet_grid(exon_biotype ~ .) +
  scale_fill_manual(values = brew) +
  scale_color_manual(values = brew) +
  guides(color = FALSE) +
  guides(fill = guide_legend("Source")) +
  xlab("Exon length (nt)")

ex_count = ggplot(data = exon, aes(x = exon_biotype, fill = paste(exon_biotype, discovery))) +
  ggtitle("Number of exons") +
  geom_bar(colour = "black",
           alpha = 0.8,
           position = "dodge") +
  scale_fill_manual(values = brew) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  xlab("Biotype") +
  ylab("Number of exons") +
  geom_text(
    stat = 'count',
    aes(label = ..count..),
    vjust = -0.2,
    position = position_dodge(width = 0.9)
  ) +
  guides(color = FALSE) +
  guides(fill = guide_legend("Source")) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_blank())

#############################################################################
# PDF
#############################################################################
pdf(paste0(prefix, ".annexa.qc.pdf"), width = 7, height = 7)
cover <- textGrob("ANNEXA report",
                  gp = gpar(fontsize = 40,
                            col = "black"))
        
sub_cover_text <- paste(
  print(Sys.Date()),
  paste("ANNEXA version:", ANNEXA_version),
  sep = "\n")

sub_cover <- textGrob(sub_cover_text,
                  gp = gpar(fontsize = 20,
                            col = "black"))

grid.arrange(cover, sub_cover)

# Command line
formatted_command <- gsub("(^|\\s)(-{1,2}\\w+)", "\n\\1\\2", command)

grid.arrange(
  textGrob("ANNEXA command line:", 
            gp = gpar(fontsize = 14, fontface = "bold")),
  textGrob(
    paste(formatted_command),
    just = "left",
    x = 0,
    hjust = 0,
    gp = gpar(fontsize = 10)
  ),
  ncol = 1,
  heights = c(0.1, 0.9)
)

# Gene
grid.arrange(textGrob(
  "Gene Characterization",
  gp = gpar(fontface = "italic", fontsize = 20),
  vjust = 0
))
count
len
iso
gene_counts
gene_validate
gene_counts_sample
gene_ext
gene_ext_dist

# Transcripts
grid.arrange(textGrob(
  "Transcript Characterization",
  gp = gpar(fontface = "italic", fontsize = 20),
  vjust = 0
))
tx_count
tx_len
tx_ex

# Exons
grid.arrange(textGrob(
  "Exon Characterization",
  gp = gpar(fontface = "italic", fontsize = 20),
  vjust = 0
))
ex_count
ex_len

dev.off()
