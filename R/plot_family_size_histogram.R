#!/usr/bin/env Rscript

# This script will plot output family size histograms as produced by the
# SWIFT UMI pipeline. It searches for files <samplename>.counts2hist and plots the
# result for all samples in one figure.

library(tidyr)
library(ggplot2) # use dev version >= 3.0.0.9000 or bear the extra Rplots.pdf
library(optparse)

option.list <- list(
  make_option(c("-i", "--input"), type = "character", default = ".",
              help = "Absolute or relative path to where the files are."),
  make_option(c("-o", "--output"), type = "character", default = "family_size_histogram",
              help = "Relative or absolute path to the output figure."),
  make_option(c("--width"), type = "double", default = 12.0,
              help = "Width of the output plot, in cm."),
  make_option(c("--height"), type = "double", default = 8.0,
              help = "Height of the output plot, in cm."),
  make_option(c("-d", "--device"), type = "character", default = "pdf",
              help = "Which device to use for plotting of the final image. Default is
                pdf, but e.g. png and jpg are also available. See ?ggsave for more
                info."))

opt <- parse_args(OptionParser(option_list = option.list))

xy <- list.files(opt$input, pattern = "^.*\\.counts2hist$", full.names = TRUE)

# Extract sample name.
samplename <- gsub("^(.*)\\.counts2hist", "\\1", basename(xy))

# Append sample to data.frame.
hs <- mapply(FUN = function(x, y) {
  d <- read.table(x)
  d$sample <- y
  colnames(d) <- c("value", "sample")
  d
}, x = xy, y = samplename, SIMPLIFY = FALSE)

# Merge all samples into one data.frame.
hs <- do.call(rbind, hs)
rownames(hs) <- NULL

# Plot and save into a figure.
out <- ggplot(hs, aes(x = value)) +
  theme_bw() +
  theme(axis.text = element_text(size = 6)) +
  ylab("Frequency") +
  xlab("MID family size") +
  geom_histogram() +
  facet_wrap(~ sample)

ggsave(filename = "Figure 1: family_size_histogram.pdf", plot = out,
       width = opt$width, height = opt$height, device = opt$dev, units = "cm",
       title = "Family size histogram")
