library(DESeq2)
library(ggplot2)
library(tidyverse)
raw_counts <- read.csv("data/raw_data/raw_counts_90min.csv", row.names = 1)
raw_counts <- as.matrix(raw_counts)
metadata <- read.csv("data/raw_data/metadata_90min.csv" , row.names = 1)

dds_matrix <- DESeqDataSetFromMatrix(countData = raw_counts,
                                     colData = metadata,
                                     design = ~treatment)
dds_matrix$treatment <- relevel(dds_matrix$treatment, ref = "control")
dds <- DESeq(dds_matrix)
rld <- rlog(dds)
plot <- plotPCA(rld, intgroup = "treatment", ntop = 500) +
  ggtitle("Treatment Conditions: PCoA") + 
  labs(color="Treatment") + 
  theme_test() +
  theme(
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key = element_rect(fill = "lightgrey"),
    panel.border = element_rect(linewidth = 2))
plot

ggsave(filename="figures/pcoa.png", width = 7, height = 5)
