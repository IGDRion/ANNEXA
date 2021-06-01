#! /usr/bin/env Rscript
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggridges)
library(RColorBrewer)

brew = c("#008387", "#a0d3db", "#ad8500", "#d9cb98")
palette = c("#00AFBB", "#FC4E07", "#E7B800")

theme_set(
  theme_pubr() +
    theme(
      panel.grid.major = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      strip.text = element_text(face = "bold")
    )
)

GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <-
      transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )
    grp <- data[1, "group"]
    newdata <-
      plyr::arrange(transform(data, x = if (grp %% 2 == 1)
        xminv
        else
          xmaxv), if (grp %% 2 == 1)
            y
        else
          - y)
    newdata <-
      rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
      round(newdata[1, "x"])
    
    if (length(draw_quantiles) > 0 &
        !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                1))
      quantiles <-
        ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <-
        data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <-
        rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <-
        GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "identity",
           ...,
           draw_quantiles = NULL,
           trim = TRUE,
           scale = "area",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomSplitViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        draw_quantiles = draw_quantiles,
        na.rm = na.rm,
        ...
      )
    )
  }


#############################################################################
# GENE
#############################################################################
gene = read.csv("gene.stats", header = T)
min_norm = min(gene[length(gene)])
gene$gene_biotype[gene$gene_biotype == "protein_coding"] = "mRNA"

# Length
med_length = gene %>%
  group_by(discovery, gene_biotype) %>%
  summarize(median = median(length), gene_biotype = gene_biotype)

len = gene %>%
  ggplot(aes(x = length, fill = paste(gene_biotype, discovery))) +
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
  guides(fill = guide_legend("Source"))


# Nb isoformes
iso = gene %>%
  mutate(isoformes = ifelse(nb_transcripts == 1, "1", "2+")) %>%
  ggplot(aes(x = discovery, fill = isoformes)) +
  geom_bar(position = "fill",
           alpha = 0.8,
           colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~ gene_biotype) +
  scale_fill_manual(values = palette)

# Nombre de gène chaque catégorie
count = ggplot(data = gene, aes(x = gene_biotype, fill = paste(gene_biotype, discovery))) +
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
  xlab("Biotype") + ylab("Number of gene")

# Quantif in ridge plot
gene_quant = gene %>%
  select(c(2, 7:dim(gene)[2])) %>%
  melt() %>%
  filter(value > min_norm) %>%
  mutate(variable = gsub("\\..*", "", variable)) %>%
  ggplot() +
  geom_density_ridges2(aes(
    x = value,
    y = variable,
    fill = paste(gene_biotype, discovery)
  ), alpha = 0.7) +
  theme(legend.position = "top") +
  facet_grid(gene_biotype ~ .) +
  scale_fill_manual(values = brew) +
  guides(fill = guide_legend("Source")) +
  xlab("Normalized Gene expression") + ylab("Sample")

# Quantif in split violin plot
gene_quant_violin = gene %>%
  select(c(2, 7:dim(gene)[2])) %>%
  melt() %>%
  filter(value > min_norm) %>%
  mutate(variable = gsub("\\..*", "", variable)) %>%
  ggplot(aes(
    x = variable,
    y = value,
    fill = paste(gene_biotype, discovery)
  )) +
  geom_split_violin(trim = T, alpha = 0.7) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = .5,
    size = .25,
    position = position_dodge(width = .5),
    show.legend = F
  ) +
  facet_grid(gene_biotype ~ .) +
  ylab("Normalized Gene Expression") +
  xlab("Sample") +
  guides(fill = guide_legend("Source")) +
  coord_flip() +
  scale_fill_manual(values = brew) +
  theme(legend.position = "top")

# 5'-3' extensions count
gene_ext = gene %>%
  filter(ext_5 > 0 | ext_3 > 0) %>%
  mutate(ext = ifelse(ext_5 > 0 & ext_3 > 0, "5'-3'",
                      ifelse(ext_5 > 0, "5'", "3'"))) %>%
  ggplot(aes(x = gene_biotype, fill = ext)) +
  geom_bar(position = "stack",
           colour = "black",
           alpha = 0.8) +
  scale_fill_manual(values = palette) +
  theme(legend.position = "top")

# 5'-3' extensions distribution
gene_ext_dist = gene %>%
  filter(ext_5 > 0 | ext_3 > 0) %>%
  select(c("ext_5", "ext_3", "gene_biotype", "discovery")) %>%
  melt() %>%
  ggplot(aes(x = value, fill = paste(gene_biotype, discovery))) +
  geom_density(alpha = 0.5) +
  facet_grid(variable ~ .) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  theme(legend.position = "top") +
  guides(fill = guide_legend("Source")) +
  xlab("Length of extension (nt)") + ylab("Density") +
  scale_fill_manual(values = c(brew[1], brew[3]))

#############################################################################
# TRANSCRIPT
#############################################################################
transcript = read.csv("transcript.stats", header = T)
transcript$transcript_biotype[transcript$transcript_biotype == "protein_coding"] = "mRNA"
lncRNA_biotypes = c("retained_intron", "lncRNA", "antisense", "non-coding")

transcript = transcript %>%
  mutate(
    transcript_biotype = if_else(
      transcript_biotype %in% lncRNA_biotypes,
      "lncRNA",
      transcript_biotype
    )
  ) %>%
  filter(transcript_biotype %in% c("mRNA", "lncRNA"))

# Length distrib
med_length_tx = transcript %>%
  group_by(discovery, transcript_biotype) %>%
  summarize(median = median(length), transcript_biotype = transcript_biotype)

tx_len = transcript %>%
  ggplot(aes(x = length, fill = paste(transcript_biotype, discovery))) +
  geom_vline(data = med_length_tx,
             aes(
               xintercept = median,
               color = paste(transcript_biotype, discovery)
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
  guides(color = FALSE)

# Mono vs multiexonique
tx_ex = transcript %>%
  mutate(exons = ifelse(nb_exons == 1, "1", "2+")) %>%
  ggplot(aes(x = discovery, fill = exons)) +
  geom_bar(position = "fill",
           alpha = 0.8,
           colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~ transcript_biotype) +
  scale_fill_manual(values = palette)

# Count
tx_count = transcript %>%
  ggplot(aes(
    x = transcript_biotype,
    fill = paste(transcript_biotype, discovery, "in", gene_discovery, "gene")
  )) +
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
  theme(legend.title = element_text(size = 10)) +
  xlab("Transcript Biotype") + ylab("Number of transcripts")

# Quantif ridge plot
tx_quant = transcript %>%
  select(c(3, 6:dim(transcript)[2])) %>%
  melt() %>%
  filter(value >= 1) %>%
  mutate(variable = gsub("\\..*", "", variable)) %>%
  ggplot() +
  geom_density_ridges2(aes(
    x = value,
    y = variable,
    fill = paste(transcript_biotype, discovery)
  ), alpha = 0.6) +
  theme(legend.position = "top") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  facet_grid(transcript_biotype ~ .) +
  scale_fill_manual(values = brew) +
  guides(fill = guide_legend("Source")) +
  xlab("Transcript expression") + ylab("Sample")

# Quantif split violin plot
tx_quant_violin = transcript %>%
  select(c(3, 6:dim(transcript)[2])) %>%
  melt() %>%
  filter(value >= 1) %>%
  mutate(variable = gsub("\\..*", "", variable)) %>%
  ggplot(aes(
    x = variable,
    y = value,
    fill = paste(transcript_biotype, discovery)
  )) +
  geom_split_violin(trim = T, alpha = 0.7) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = .5,
    size = .25,
    position = position_dodge(width = .5),
    show.legend = F
  ) +
  facet_grid(transcript_biotype ~ .) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                 .x))
  ) +
  ylab("Transcript Expression") +
  xlab("Sample") +
  coord_flip() +
  theme(legend.position = "top") +
  guides(fill = guide_legend("Source")) +
  scale_fill_manual(values = brew)

#############################################################################
# EXON
#############################################################################
exon = read.csv("exon.stats", header = T)
exon$exon_biotype[exon$exon_biotype == "protein_coding"] = "mRNA"

# Length distrib
med_length_ex = exon %>%
  group_by(discovery, exon_biotype) %>%
  summarize(median = median(length), exon_biotype = exon_biotype)

ex_len = exon %>%
  ggplot(aes(x = length, fill = paste(exon_biotype, discovery))) +
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
  guides(fill = guide_legend("Source"))

ex_count = ggplot(data = exon, aes(x = exon_biotype, fill = paste(exon_biotype, discovery))) +
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
  theme(legend.position = "top")

#############################################################################
# PDF
#############################################################################
pdf("qc_gtf.pdf", width = 11, height = 8)
ggarrange(
  len,
  ggarrange(
    iso,
    count,
    ncol = 2,
    labels = c("B", "C"),
    legend = "left",
    common.legend = TRUE
  ),
  nrow = 2,
  labels = "A",
  legend = "top"      # Étiquette du line plot
) %>%
  annotate_figure(top = text_grob("Gene",
                                  face = "bold",
                                  size = 18))
ggarrange(
  tx_len,
  ggarrange(
    tx_ex,
    tx_count,
    ncol = 2,
    labels = c("B", "C"),
    legend = "left",
    common.legend = F
  ),
  nrow = 2,
  labels = "A",
  legend = "top"      # Étiquette du line plot
) %>%
  annotate_figure(top = text_grob("Transcript",
                                  face = "bold",
                                  size = 18))

ggarrange(
  ex_len,
  ggarrange(
    ex_count,
    ggplot() + theme_void(),
    ncol = 2,
    labels = "B",
    legend = "left",
    common.legend = F
  ),
  nrow = 2,
  labels = "A",
  legend = "top"      # Étiquette du line plot
) %>%
  annotate_figure(top = text_grob("Exon",
                                  face = "bold",
                                  size = 18))

gene_quant %>%
  annotate_figure(top = text_grob("Normalized Gene expression",
                                  face = "bold",
                                  size = 18))
gene_quant_violin %>%
  annotate_figure(top = text_grob("Normalized Gene expression",
                                  face = "bold",
                                  size = 18))
tx_quant %>%
  annotate_figure(top = text_grob("Transcript expression",
                                  face = "bold",
                                  size = 18))
tx_quant_violin %>%
  annotate_figure(top = text_grob("Transcript expression",
                                  face = "bold",
                                  size = 18))
ggarrange(
  gene_ext_dist,
  ggarrange(
    gene_ext,
    ggplot() + theme_void(),
    ncol = 2,
    labels = "B",
    legend = "left",
    common.legend = F
  ),
  nrow = 2,
  labels = "A",
  legend = "top"
) %>%
  annotate_figure(top = text_grob("Gene extensions",
                                  face = "bold",
                                  size = 18))

dev.off()
