library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(rstudioapi)

# Get the directory of the current script
script_path <- normalizePath(rstudioapi::getSourceEditorContext()$path)
script_dir <- dirname(script_path)

# Set the working directory
setwd(script_dir)

# Input file
dat <- read.csv("./output_canFam4_ref_final.csv", header = TRUE)

# Reformat input with feat_level both for gn and tx
dat <- dat %>%
  pivot_longer(cols = c(Gene, Tx),
               names_to = "feat_level", values_to = "count") %>%
  mutate(
    feat_level = replace(feat_level, feat_level == "Tx", "Transcripts"),
    Full_vs_Filter = case_when(
      Full_vs_Filter == "Full" ~ "Reference Annotation",
      Full_vs_Filter == "Filter" ~ "Reconstructed Annotation",
      TRUE ~ Full_vs_Filter
    ))

# Set the order of Tool levels
dat$Tool <- factor(dat$Tool, levels = c("Ref", "Bambu", "Stringtie"))

p <- dat %>%
  ggplot(aes(x = Tool, y = count,
             fill =  Annotation, alpha = Full_vs_Filter)) +
  geom_col(position = "identity") +
  geom_text(aes(label = count, alpha = NULL), 
            position = position_stack(vjust = 1.05),
            size = 3.5) +
  scale_fill_brewer(palette = "Set1") +
  scale_alpha_discrete(range = c(0.7, 0.25),
                       labels = c("Reconstructed", "Reference (Ref)"),
                       name = "Annotation Type") +
  scale_color_manual(values = c("black")) +
  facet_grid(feat_level ~ Annotation,  scales = "free_y") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  scale_x_discrete(limits = c("Ref", "Bambu", "Stringtie")) +
  theme_minimal(base_size = 16) +
  xlab("canFam4") +
  ylab("Number of elements (in thousands)") +
  guides(fill = "none", color = "none",
         alpha = guide_legend(override.aes = list(label = ""))) +
  theme(
    # remove the vertical grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(linewidth = 0.1, color = "black"),
    # Size of font
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11)
  )

# Save as PDF
p
pdf("canFam4_ref.pdf", width = 10, height = 8)
print(p)
dev.off()