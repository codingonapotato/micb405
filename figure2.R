# Define locus tags for genes of interest
efflux <- c("HKO16_RS09100", "HKO16_RS09095", "HKO16_RS12235", "HKO16_RS12240", 
           "HKO16_RS12245", "HKO16_RS14440", "HKO16_RS14445", "HKO16_RS14450", 
           "HKO16_RS12230", "HKO16_RS10335", "HKO16_RS09105", "HKO16_RS09110", 
           "HKO16_RS09120", "HKO16_RS16260")

lactamase <- c("HKO16_RS12545", "HKO16_RS07685")

genes <- c(efflux, lactamase)

# Read in DESeq2 results for carbapenem treatments into proper format for plotting
master <- data.frame()
for (file in list.files("results/filtered", full.names = TRUE, pattern = ".*penem.*")) {
  cname <- str_remove_all(file, "(.*/|_sig|\\.tsv)")
  
  df <- read.csv(file, sep = "\t") %>%
    select(tag, "{cname}" := log2FoldChange) %>%
    dplyr::filter(tag %in% genes)
  
  if (nrow(master) == 0) {
    master <- df
  } else {
    master <- full_join(master, df)
  }
}

# Read in gene annotations and join with DESeq2 results
meta <- read.csv("data/genes_simple.tsv", sep = '\t') %>% select(-type, -protein)
master <- left_join(master, meta)


# Additional data formatting for plotting
imi_data <- select(master, tag, symbol,
                   "30min" = imipenem_30,
                   "90min" = imipenem_90) %>% 
  arrange(desc(`90min`))
imi_labels <- paste(imi_data$symbol, imi_data$tag, sep = "  ")
imi_data <- select(imi_data, -tag, -symbol) %>% as.matrix()

mero_data <- select(master, tag, symbol,
                   "30min" = meropenem_30,
                   "90min" = meropenem_90) %>% 
  arrange(desc(`90min`))
mero_labels <- paste(mero_data$symbol, mero_data$tag, sep = "  ")
mero_data <- select(mero_data, -tag, -symbol) %>% as.matrix()

# Define colors for heatmap
colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
na_color <- "grey"

# Plot heatmaps
imi_hm <- Heatmap(
            imi_data,
            name = "Log2FC",
            na_col = na_color,           
            col = colors,           
            cluster_rows = FALSE,        
            cluster_columns = FALSE,     
            row_names_side = "left",     
            column_names_side = "top",   
            row_labels = imi_labels,
            column_names_rot = 360,
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            width = unit(3, "cm"),
            height = unit(7, "cm")
          )

mero_hm <- Heatmap(
              mero_data,
              name = "Log2FC",
              na_col = na_color,         
              col = colors,         
              cluster_rows = FALSE,      
              cluster_columns = FALSE,   
              row_names_side = "left",   
              column_names_side = "top", 
              row_labels = mero_labels,
              column_names_rot = 360,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              width = unit(3, "cm"),
              height = unit(7, "cm")
            )

imi_hm
mero_hm


# Plot volcano plots
dat <- read.csv("results/90min/imipenem_vs_ctrl_90.csv")
dat_annotated <- dat %>%
  mutate(Annotation  = case_when(
    gene_id %in% genes & padj < 0.05 ~ "GOI",
    TRUE ~ "Other"
  )) %>%
  arrange(Annotation == "GOI")
dat_annotated %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = Annotation)) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("darkorange", "azure2")) 

