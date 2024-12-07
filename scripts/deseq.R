# Read in raw counts files and convert to matrix
raw_counts <- read.csv("data/raw_counts_90min.csv", row.names = 1)
raw_counts <- as.matrix(raw_counts)

# Read in metadata file
metadata <- read.csv("data/metadata_90min.csv" , row.names = 1)

# Create DESeq2 obj and set control condition
dds_matrix <- DESeqDataSetFromMatrix(countData = raw_counts,
                                     colData = metadata,
                                     design = ~treatment)
dds_matrix$treatment <- relevel(dds_matrix$treatment, ref = "control")

# Run DESeq2
dds <- DESeq(dds_matrix)

# Define antibiotics and read in gene annotations
antibiotics <- c("ciprofloxacin", "colistin", "imipenem", "meropenem", "moxifloxacin", "tobramycin")
genes <- read.csv("data/genes_simple.tsv", sep = "\t")

# Output DESeq2 results for each treatment-control comparison (filtered and annotated)
for (ab in antibiotics) {
  res <- results(dds, name = glue("treatment_{ab}_vs_control")) %>% as.data.frame()
  filtered <- dplyr::filter(res, padj <= 0.05)
  filtered$tag <- row.names(filtered)
  row.names(filtered) <- NULL
  joined <- dplyr::left_join(filtered, genes, by = 'tag')
  write.table(joined, glue("results/filtered/{ab}_sig_90.tsv"), quote = FALSE, row.names = FALSE, sep = '\t')
}

# Output raw DESeq2 results for each comparison
for (ab in antibiotics) {
  res <- results(dds, name = glue("treatment_{ab}_vs_control")) %>% as.data.frame()
  temp <- cbind(rownames(res), data.frame(res, row.names=NULL)) %>% dplyr::rename(gene_id = 1)
  write.csv(temp, glue("results/30min/{ab}_vs_ctrl_90.csv"), row.names = FALSE, quote = FALSE)
}




