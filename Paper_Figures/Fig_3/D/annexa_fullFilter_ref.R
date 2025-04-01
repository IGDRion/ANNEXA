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


# INput file
dat <- read.csv("./output.csv", header = TRUE)

# Reformet input with feat_level both for gn and tx
dat <- dat %>%
  pivot_longer(cols = c(Gene, Tx),
               names_to = "feat_level", values_to = "count") %>%
  mutate(
        feat_level = replace(feat_level, feat_level == "Tx", "Transcripts"))

p <- dat %>%
  ggplot(aes(x = Tool, y = count,
             fill =  Annotation, alpha = Full_vs_Filter)) +
  geom_col(position = "identity") +
  scale_fill_manual(values = c("#964B00", "#7A288A")) +
  scale_alpha_discrete(range = c(0.7, 0.25),
                       labels = c("Filter", "Full"),
                       name = "Filtering Operation") +
  facet_grid(feat_level ~ Annotation,  scales = "free_y") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  theme_minimal(base_size = 16) +
  xlab("GRCh38") +
  ylab("Number of elements (in thousands)") +
  guides(fill = "none") +
  theme(
    # remove the vertical grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(linewidth = 0.1, color = "black"),
    # Size of font
    strip.text.x = element_text(size = 16, face = "bold"),
    strip.text.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16)
  )

p
pdf("human.pdf")
print(p)
dev.off()