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
dat <- read.csv("./biotype_full_info.csv", header = TRUE)

# Reformat input with feat_level both for gn and tx
dat <- dat %>%
  pivot_longer(cols = c(Gene_coding,Gene_lncRNA,TX_Coding,TX_lncRNA),
               names_to = "feat_level", values_to = "count") %>%
  mutate(
    feat_level = case_when(
      feat_level == "Gene_coding" ~ "mRNA Gene",
      feat_level == "Gene_lncRNA" ~ "lncRNA Gene",
      feat_level == "TX_Coding" ~ "mRNA Tx",
      feat_level == "TX_lncRNA" ~ "lncRNA Tx",
      TRUE ~ feat_level))

# Reorder feat_level to ensure Coding is on the left and Non Coding is on the right
dat$feat_level <- factor(dat$feat_level, 
                         levels = c("mRNA Gene", "lncRNA Gene", 
                                    "mRNA Tx", "lncRNA Tx"))

# Reorder Tool to ensure columns are in Ref, Bambu, Stringtie order
dat$Tool <- factor(dat$Tool, levels = c("Ref", "Bambu", "Stringtie"))

# Reshape data for overlaid plot
dat_reshaped <- dat %>%
  pivot_wider(names_from = Tool, values_from = count) %>%
  pivot_longer(cols = c(Bambu, Stringtie), 
               names_to = "Tool", 
               values_to = "count")

dat_reshaped$Tool <- factor(dat_reshaped$Tool, levels = c("Bambu", "Stringtie"))

# Create the plot
p <- ggplot(dat_reshaped, aes(x = Tool, y = count)) +
  geom_col(aes(y = Ref, alpha = "Reference"), fill = "#964B00") +
  geom_col(aes(alpha = "Reconstructed"), fill = "#964B00") +
  scale_alpha_manual(values = c("Reference" = 0.25, "Reconstructed" = 0.7),
                     name = "Annotation Type") +
  facet_wrap(~ feat_level, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"),
                     limits = c(0, NA)) +
  theme_minimal(base_size = 16) +
  xlab("GRCh38 - Gencode") +
  ylab("Number of elements (in thousands)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.1, color = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_cartesian(ylim = c(0, NA))

# Display the plot
p

# Save the plot
pdf("human_overlaid.pdf")
print(p)
dev.off()
