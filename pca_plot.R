library(DESeq2)
library(ggplot2)
library(tidyverse)
raw_counts <- read.csv("data/raw_counts_90min.csv", row.names = 1)
raw_counts <- as.matrix(raw_counts)
metadata <- read.csv("data/metadata_90min.csv" , row.names = 1)
antibiotics <- c("ciprofloxacin", "colistin", "imipenem", "meropenem", "moxifloxacin", "tobramycin")
targets <- c("DNA gyrase", "Membrane", "Cell wall", "Cell wall", "DNA gyrase", "Ribosome")
ab_target_map <- data.frame(antibiotics, targets)
augmented_metadata <- join()

dds_matrix <- DESeqDataSetFromMatrix(countData = raw_counts,
                                     colData = metadata,
                                     design = ~treatment)
dds_matrix$treatment <- relevel(dds_matrix$treatment, ref = "control")
dds <- DESeq(dds_matrix)
rld <- rlog(dds)
plot <- plotPCA(rld, intgroup = "treatment", ntop = 500) +
  ggtitle("Treatment Conditions: PCoA") + 
  labs(color="Treatments")
plot
