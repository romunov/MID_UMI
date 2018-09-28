#!/usr/bin/env Rscript

# This script will plot output from SWIFT UMI pipeline.
# Input file is variant_AFbins_prePOST-MID.txt. Run the script using e.g.
# Rscript ./R/plot_variant_af_bins.R --input variant_AFbins_prePOST-MID.txt

library(tidyr)
library(ggplot2) # use dev version >= 3.0.0.9000 or bear the extra Rplots.pdf
library(optparse)

option.list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file, normally variant_AFbins_prePOST-MID from
                the pipeline"),
  make_option(c("-o", "--output"), type = "character", default = "variant_AFbins_prePOST-MID",
              help = "Relative or absolute path to the output figure."),
  make_option(c("--width"), type = "double", default = 10.0,
              help = "Width of the output plot, in cm."),
  make_option(c("--height"), type = "double", default = 10.0,
              help = "Height of the output plot, in cm."),
  make_option(c("-d", "--device"), type = "character", default = "pdf",
              help = "Which device to use for plotting of the final image. Default is
                pdf, but e.g. png and jpg are also available. See ?ggsave for more
                info."))

opt <- parse_args(OptionParser(option_list = option.list))

# Import data and reflow it into a long format which is handy for plotting.
xy <- read.table(opt$input, header = TRUE, sep = "\t")
xy <- droplevels(xy[-1, ]) # remove first line because it's a summary

# Save this to specify the correct order of levels.
lvl.order <- as.character(xy$Sample)

# Reflow data into long format.
xy <- gather(xy, key = tmp, value = value, -Sample)
colnames(xy) <- c("vars", "tmp", "value")
xy$vars <- factor(xy$vars, levels = lvl.order, ordered = TRUE)

# Extract sample and origin (fgbio or preMID) information and store it in a new variable.
xy$sample <- gsub("^(.*)(\\.fgbio|\\.preMID$)", "\\1", xy$tmp)
# fgbio is the post MID treatment.
xy$origin <- factor(gsub("^(.*\\.)(fgbio|preMID$)", "\\2", xy$tmp), 
                    levels = c("fgbio", "preMID"),
                    labels = c("postMID", "preMID"))
xy <- xy[, c("vars", "sample", "origin", "value")]

# Calculate summary stats per sample and origin.
xy.stats <- aggregate(value ~ origin + sample, FUN = sum, data = xy)

out <- ggplot(xy, aes(x = vars, y = value, fill = origin)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text = element_text(size = 6),
        legend.position = "top") +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  geom_col(position = "dodge", width = 0.75) +
  xlab("Variant allele frequencies") + ylab("Count") +
  facet_wrap(~ sample, ncol = 5)

ggsave(paste("Figure2", paste(opt$output, opt$dev, sep = "."), sep = "_"),
       width = opt$width, height = opt$height, device = opt$dev, units = "cm",
       title = "Variant allele frequencies")
