library(tidyr)
library(ggplot2)

# Import data and reflow it into a long format which is handy for plotting.
xy <- read.table("../variant_AFbins_prePOST-MID.txt", header = TRUE, sep = "\t")
xy <- droplevels(xy[-1, ]) # remove first line because it's a summary

# Save this to specify the correct order of levels.
lvl.order <- as.character(xy$Sample)

xy <- gather(xy, key = tmp, value = value, -Sample)
colnames(xy) <- c("vars", "tmp", "value")
xy$vars <- factor(xy$vars, levels = lvl.order, ordered = TRUE)

# Extract sample and origin (fgbio or preMID) information and store it in a new variable.
xy$sample <- gsub("^(.*)(\\.fgbio|\\.preMID$)", "\\1", xy$tmp)
xy$origin <- gsub("^(.*\\.)(fgbio|preMID$)", "\\2", xy$tmp)
xy <- xy[, c("vars", "sample", "origin", "value")]

# Calculate summary stats per sample and origin.
xy.stats <- aggregate(value ~ origin + sample, FUN = sum, data = xy)

ggplot(xy, aes(x = vars, y = value, fill = origin)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "top") +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  geom_col(position = "dodge", width = 0.75) +
  xlab("Vars AF") + ylab("Count") +
  facet_wrap(~ sample, ncol = 5)

ggsave("../variant_AFbins_prePOST-MID.pdf", width = 6, height = 6)
