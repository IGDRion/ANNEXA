library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(rstudioapi)

script_path <- normalizePath(rstudioapi::getSourceEditorContext()$path)
script_dir <- dirname(script_path)
setwd(script_dir)

x <- read.csv2("output_canFam4.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

#x <- x[grep("ensembl", x$file, ignore.case = TRUE), ]

# Separate file column into 4 columns
x_sep <- x %>% separate(file,
                        into = c("annot", "type", "tool", "gffcmp_lvl"),
                        sep = "_",
                        convert = TRUE, extra = "merge") %>%
  mutate(gffcmp_lvl = factor(gffcmp_lvl, levels = c("exon", "transcript", "locus")))

# pivot longer metrics
x_sep_format <- x_sep %>% 
  pivot_longer(c("Sensitivity", "Precision"),
               names_to = "perf_metrics",
               values_to = "perf_values") %>%
  mutate_at("perf_values", as.numeric)


# create Figure
figure <- ggplot(data = x_sep_format,
                 aes(x = perf_values,
                     y = paste("ANNEXA", type),
                     color = perf_metrics,
                     shape = as.factor(annot))) +
  geom_point(size = 6) +
  labs(
       title = "Benchmark Stringtie vs Bambu wrt Annotation - canFam 4",
       y = " ",
       x = "Performance") +
  scale_color_discrete(name = "Performance Metrics") +
  scale_shape_discrete(name = "Annotation") +
  theme_bw(base_size = 18) +
  facet_grid(tool ~ gffcmp_lvl ) +
  theme(strip.text.y = element_text(angle = 0))


plot(figure)

pdf("plot_benchmark.pdf", width = 14)
plot(figure)
dev.off()
