# Create a master raw counts file from STAR gene counts files.
aligned <- "PRJNA686448/aligned2"

master <- data.frame(dummy = c(1:3777))

for (dir in list.dirs(aligned, recursive = FALSE, full.names = FALSE)) {
  geneCounts <- read.csv(paste0(aligned, "/", dir, "/", "ReadsPerGene.out.tab"), sep = '\t', skip = 3)
  vec <- geneCounts[, 4]
  df <- data.frame(vec)
  colnames(df) <- dir
  master <- cbind(master, df)
}

master$dummy <- geneCounts[,1]
colnames(master)[colnames(master) == 'dummy'] <- 'gene'


# Create master meta data file with proper formatting
metadata <- "PRJNA686448/sra_explorer_metadata.tsv"
meta <- read.csv(metadata, sep="\t")

meta <- mutate(meta,
               treatment = case_when(grepl("control", Title, fixed = TRUE) ~ "control",
                                     TRUE ~ str_split_i(Title, " ", 5)),
               timepoint = str_split_i(Title, " ", -1)) %>% 
  select(treatment, timepoint, Accession, Title) %>%
  mutate(sample = paste(treatment, timepoint, sep = "_")) %>%
  arrange(sample)

# Further formatting of metadata for compatibility with DESeq2 analysis
curr <- ""
count <- 1
for(i in 1:nrow(meta)) {
  if (meta[i,5] != curr) {
    curr <- meta[i,5]
    count <- 1
    meta[i,5] = paste(meta[i,5], count, sep = "_")
  } else {
    meta[i,5] <- paste(meta[i,5], count, sep = "_")
  }
  count <- count + 1
}
row.names(meta) <- meta$sample

# Divide meta data by time point
meta_90min <- filter(meta, timepoint == "90min") %>% arrange(Accession) %>% select(treatment)
meta_30min <- filter(meta, timepoint == "30min") %>% arrange(Accession) %>% select(treatment)

write.csv(meta_30min, "data/metadata_30min.csv", quote = FALSE)
write.csv(meta_90min, "data/metadata_90min.csv", quote = FALSE)

# Formatting raw counts file for compatibility with DESeq2 + separate by timepoint
cols <- c("gene", row.names(arrange(meta, Accession)))
colnames(master) <- cols

meta <- arrange(meta, Accession) %>% select(treatment, timepoint)
write.csv(meta, "data/metadata.csv", quote = FALSE)

master_90min <- select(master, contains("90min"))
row.names(master_90min) <- master$gene
master_30min <- select(master, contains("30min"))
row.names(master_30min) <- master$gene

row.names(master) <- master$gene
master <- select(master, -gene)

write.csv(master, "data/raw_counts.csv", quote = FALSE)
write.csv(master_30min, "data/raw_counts_30min.csv", quote = FALSE)
write.csv(master_90min, "data/raw_counts_90min.csv", quote = FALSE)


# Create a simplified gene annotation file
genes <- read.csv("data/genes.tsv", sep = '\t')
genes <- select(genes, symbol = Symbol, type = Gene.Type, protein = Name, tag = Locus.tag)

# Create raw counts and metadata files for time point comparisons
raw_counts <- read.csv("data/raw_counts.csv", row.names = 1)
control_counts <- select(raw_counts, contains("control"))
mero_counts <- select(raw_counts, contains("meropenem"))

meta <- read.csv("data/metadata.csv", row.names = 1)
mero_meta <- dplyr::filter(meta, treatment == "meropenem") %>% select(timepoint)
ctrl_meta <- dplyr::filter(meta, treatment == "control") %>% select(timepoint)

write.csv(control_counts, "time_data/control_counts.csv", quote = FALSE, row.names = TRUE)
write.csv(mero_counts, "time_data/mero_counts.csv", quote = FALSE, row.names = TRUE)
write.csv(ctrl_meta, "time_data/control_meta.csv", quote = FALSE, row.names = TRUE)
write.csv(mero_meta, "time_data/mero_meta.csv", quote = FALSE, row.names = TRUE)


# Data for time comparisions by separate by treatment
treatment <- c("ciprofloxacin", "colistin", "imipenem", "meropenem", "moxifloxacin", "tobramycin", "control")
meta <- read.csv("data/metadata.csv", row.names = 1)
counts <- read.csv("data/raw_counts.csv", row.names = 1)

for (t in treatment) {
  meta_df <- dplyr::filter(meta, treatment == t) %>% select(timepoint)
  counts_df <- dplyr::select(counts, contains(t))
  write.csv(meta_df, glue("time_data/all/{t}_meta.csv"), quote = FALSE, row.names = TRUE)
  write.csv(counts_df, glue("time_data/all/{t}_counts.csv"), quote = FALSE, row.names = TRUE)
}


