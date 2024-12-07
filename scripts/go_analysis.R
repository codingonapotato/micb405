suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))

library(stringr) 
library(GO.db)

#### Load GO mapping file and create character vector with all gene IDs from the file ####
geneID2GO <- readMappings("go_mapping.tsv")
geneUniverse <- names(geneID2GO)

ensure_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(sprintf("Directory created: %s", dir))
  } else {
    message(sprintf("Directory exists: %s", dir))
  }
}

# Define enums
antibiotics <- list(
  CIPROFLOXACIN = "ciprofloxacin",
  MOXIFLOXACIN = "moxifloxacin",
  TOBRAMYCIN = "tobramycin",
  IMIPENEM = "imipenem",
  MEROPENEM = "meropenem",
  COLISTIN = "colistin",
  CONTROL = "control"
)

experiments <- list(
  VS_CTRL_30 = "vs_ctrl_30",
  VS_CTRL_90 = "vs_ctrl_90",
  VS_30_90 = "30_vs_90"
)

times <- list(
  T30M = "30m",
  T90M = "90m"
)

antibiotic1 <- antibiotics$TOBRAMYCIN
antibiotic2 <- antibiotics$CIPROFLOXACIN
antibiotic3 <- antibiotics$MOXIFLOXACIN

exp <- experiments$VS_CTRL_90
time <- times$T90M

#### Load differential expression analysis data ####
csvDir = "csv"
ensure_dir(csvDir)
plotDir = "plots"
ensure_dir(plotDir)
file1 = paste(csvDir, paste(paste(antibiotic1, exp, sep="_"), "csv", sep="."), sep="/")
file2 = paste(csvDir, paste(paste(antibiotic2, exp, sep="_"), "csv", sep="."), sep="/")
file3 = paste(csvDir, paste(paste(antibiotic3, exp, sep="_"), "csv", sep="."), sep="/")
colistin_vs_control_90 <- read_csv(file1)
imipenem_vs_control_90 <- read_csv(file2)
meropenem_vs_control_90 <- read_csv(file3)

#### Find statistically significant up/down-regulation of genes ####
# Upregulated genes
colistin_up_genes <- colistin_vs_control_90 %>%
  filter(padj <= 0.05 & log2FoldChange >= 0)
imipenem_up_genes <- imipenem_vs_control_90 %>%
  filter(padj <= 0.05 & log2FoldChange >= 0)
meropenem_up_genes <- meropenem_vs_control_90 %>%
  filter(padj <= 0.05 & log2FoldChange >= 0)
# Downregulated genes
colistin_down_genes <- colistin_vs_control_90 %>%
  filter(padj <= 0.05 & log2FoldChange <= 0)
imipenem_down_genes <- imipenem_vs_control_90 %>%
  filter(padj <= 0.05 & log2FoldChange <= 0)
meropenem_down_genes <- meropenem_vs_control_90 %>%
  filter(padj <= 0.05 & log2FoldChange <= 0)

#### Extract up/down-regulated genes for each drug condition into a character vector ####
colistin_upregulated_genes <- colistin_up_genes$gene_id %>% as.character()
imipenem_upregulated_genes <- imipenem_up_genes$gene_id %>% as.character()
meropenem_upregulated_genes <- meropenem_up_genes$gene_id %>% as.character()

colistin_downregulated_genes <- colistin_down_genes$gene_id %>% as.character()
imipenem_downregulated_genes <- imipenem_down_genes$gene_id %>% as.character()
meropenem_downregulated_genes <- meropenem_down_genes$gene_id %>% as.character()

#### Generate vector of integers marking if a gene from the mapping file is present (1) or not (0) in vector of up/down-regulated genes ####
colistin_up_gene_list <- factor(as.integer(geneUniverse %in% colistin_upregulated_genes))
imipenem_up_gene_list <- factor(as.integer(geneUniverse %in% imipenem_upregulated_genes))
meropenem_up_gene_list <- factor(as.integer(geneUniverse %in% meropenem_upregulated_genes))

colistin_down_gene_list <- factor(as.integer(geneUniverse %in% colistin_downregulated_genes))
imipenem_down_gene_list <- factor(as.integer(geneUniverse %in% imipenem_downregulated_genes))
meropenem_down_gene_list <- factor(as.integer(geneUniverse %in% meropenem_downregulated_genes))

#### Set names for gene list ####
names(colistin_up_gene_list) <- geneUniverse
names(imipenem_up_gene_list) <- geneUniverse
names(meropenem_up_gene_list) <- geneUniverse

names(colistin_down_gene_list) <- geneUniverse
names(imipenem_down_gene_list) <- geneUniverse
names(meropenem_down_gene_list) <- geneUniverse

abTime1 = paste(str_to_title(antibiotic1), time, sep="_")
abTime2 = paste(str_to_title(antibiotic2), time, sep="_")
abTime3 = paste(str_to_title(antibiotic3), time, sep="_")

#### Build GOdata objects for up/down-regulated data ####
colistin_up_GO_data <- new("topGOdata", 
                  description = abTime1, 
                  ontology = "BP", 
                  allGenes = colistin_up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)
imipenem_up_GO_data <- new("topGOdata", 
                           description = abTime2, 
                           ontology = "BP", 
                           allGenes = imipenem_up_gene_list,
                           annot = annFUN.gene2GO,
                           gene2GO = geneID2GO)
meropenem_up_GO_data <- new("topGOdata", 
                           description = abTime3, 
                           ontology = "BP", 
                           allGenes = meropenem_up_gene_list,
                           annot = annFUN.gene2GO,
                           gene2GO = geneID2GO)

colistin_down_GO_data <- new("topGOdata",
                    description = abTime1,
                    ontology = "BP",
                    allGenes = colistin_down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
imipenem_down_GO_data <- new("topGOdata",
                             description = abTime2,
                             ontology = "BP",
                             allGenes = imipenem_down_gene_list,
                             annot = annFUN.gene2GO,
                             gene2GO = geneID2GO)
meropenem_down_GO_data <- new("topGOdata",
                             description = abTime3,
                             ontology = "BP",
                             allGenes = meropenem_down_gene_list,
                             annot = annFUN.gene2GO,
                             gene2GO = geneID2GO)

#### Perform Fisher's exact test with weight01 algorithm to detect enriched genes ####
colistin_up_result <- runTest(colistin_up_GO_data,
                              algorithm = "weight01",
                              statistic = "fisher")
imipenem_up_result <- runTest(imipenem_up_GO_data,
                              algorithm = "weight01",
                              statistic = "fisher")
meropenem_up_result <- runTest(meropenem_up_GO_data,
                              algorithm = "weight01",
                              statistic = "fisher")

colistin_down_result <- runTest(colistin_down_GO_data,
                                algorithm = "weight01",
                                statistic = "fisher")
imipenem_down_result <- runTest(imipenem_down_GO_data,
                                algorithm = "weight01",
                                statistic = "fisher")
meropenem_down_result <- runTest(meropenem_down_GO_data,
                                algorithm = "weight01",
                                statistic = "fisher")

upRes1 = paste(antibiotic1, "up_result", sep="_")
upRes2 = paste(antibiotic2, "up_result", sep="_")
upRes3 = paste(antibiotic3, "up_result", sep="_")

downRes1 = paste(antibiotic1, "down_result", sep="_")
downRes2 = paste(antibiotic2, "down_result", sep="_")
downRes3 = paste(antibiotic3, "down_result", sep="_")

#### Extract a summary of up/down-regulated gene results ####
colistin_up_GO <- GenTable(colistin_up_GO_data,
                           weight01 = colistin_up_result,
                           orderBy = upRes1,
                           ranksOf = upRes1,
                           topNodes = 50)
imipenem_up_GO <- GenTable(imipenem_up_GO_data,
                           weight01 = imipenem_up_result,
                           orderBy = upRes2,
                           ranksOf = upRes2,
                           topNodes = 50)
meropenem_up_GO <- GenTable(meropenem_up_GO_data,
                           weight01 = meropenem_up_result,
                           orderBy = upRes3,
                           ranksOf = upRes3,
                           topNodes = 50)

colistin_down_GO <- GenTable(colistin_down_GO_data,
                             weight01 = colistin_down_result,
                             orderBy = downRes1,
                             ranksOf = downRes1,
                             topNodes = 50)
imipenem_down_GO <- GenTable(imipenem_down_GO_data,
                             weight01 = imipenem_down_result,
                             orderBy = downRes2,
                             ranksOf = downRes2,
                             topNodes = 50)
meropenem_down_GO <- GenTable(meropenem_down_GO_data,
                             weight01 = meropenem_down_result,
                             orderBy = downRes3,
                             ranksOf = downRes1,
                             topNodes = 50)

#### Filter out any non-significant data and calculate the gene ratio ####
colistin_up_GO_filtered <- colistin_up_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)
imipenem_up_GO_filtered <- imipenem_up_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)
meropenem_up_GO_filtered <- meropenem_up_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

colistin_down_GO_filtered <- colistin_down_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)
imipenem_down_GO_filtered <- imipenem_down_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)
meropenem_down_GO_filtered <- meropenem_down_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

#### Generate visualizations ####
# Arrange the data based on the enrichment ratio. 
colistin_up_GO_filtered_arranged <- colistin_up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))
imipenem_up_GO_filtered_arranged <- imipenem_up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))
meropenem_up_GO_filtered_arranged <- meropenem_up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

colistin_down_GO_filtered_arranged <- colistin_down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))
imipenem_down_GO_filtered_arranged <- imipenem_down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))
meropenem_down_GO_filtered_arranged <- meropenem_down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Extract the order of the term column
colistin_up_order_term <- colistin_up_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector
imipenem_up_order_term <- imipenem_up_GO_filtered_arranged %>% 
  pull(Term) 
meropenem_up_order_term <- meropenem_up_GO_filtered_arranged %>% 
  pull(Term)

colistin_down_order_term <- colistin_down_GO_filtered_arranged %>% 
  pull(Term)
imipenem_down_order_term <- imipenem_down_GO_filtered_arranged %>% 
  pull(Term) 
meropenem_down_order_term <- meropenem_down_GO_filtered_arranged %>% 
  pull(Term)

# Generate ggplots
gg_colistin_up <- colistin_up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = colistin_up_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # changes the scale of the axes
gg_imipenem_up <- imipenem_up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = imipenem_up_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # changes the scale of the axes
gg_meropenem_up <- meropenem_up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = meropenem_up_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # changes the scale of the axes

gg_colistin_down <- colistin_down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = colistin_down_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # changes the scale of the axes
gg_imipenem_down <- imipenem_down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = imipenem_down_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # changes the scale of the axes
gg_meropenem_down <- meropenem_down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = meropenem_down_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # changes the scale of the axes


outUp1 = paste(plotDir, paste(paste(upRes1, exp, sep="_"), "png", sep="."), sep="/")
outUp2 = paste(plotDir, paste(paste(upRes2, exp, sep="_"), "png", sep="."), sep="/")
outUp3 = paste(plotDir, paste(paste(upRes3, exp, sep="_"), "png", sep="."), sep="/")

outDown1 = paste(plotDir, paste(paste(downRes1, exp, sep="_"), "png", sep="."), sep="/")
outDown2 = paste(plotDir, paste(paste(downRes2, exp, sep="_"), "png", sep="."), sep="/")
outDown3 = paste(plotDir, paste(paste(downRes3, exp, sep="_"), "png", sep="."), sep="/")

#### Save visualizations to file ####
ggsave(filename=outUp1, plot = gg_colistin_up,
       width = 15, height = 7)
ggsave(filename=outUp2, plot = gg_imipenem_up,
       width = 15, height = 7)
ggsave(filename=outUp3, plot = gg_meropenem_up,
       width = 15, height = 7)

ggsave(filename=outDown1, plot = gg_colistin_down,
       width = 15, height = 7)
ggsave(filename=outDown2, plot = gg_imipenem_down,
       width = 15, height = 7)
ggsave(filename=outDown3, plot = gg_meropenem_down,
       width = 15, height = 7)



### Shared Pathways ###

# Read the GTF file
gtf_file <- "GCF_013372085.1_ASM1337208v1_genomic.gtf"

# Load the GTF file as a data frame
gtf_data <- read.table(
  gtf_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

# Add column names
colnames(gtf_data) <- c(
  "seqname", "source", "feature", "start", "end", "score", 
  "strand", "frame", "attributes"
)

# Extract gene-protein mappings
extract_attribute <- function(attr_string, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  match <- regmatches(attr_string, regexec(pattern, attr_string))
  sapply(match, function(x) ifelse(length(x) > 1, x[2], NA))
}


gtf_data$gene_id <- extract_attribute(gtf_data$attributes, "gene_id")
gtf_data$product <- extract_attribute(gtf_data$attributes, "product")

# Filter rows with both gene_id and product
gene_protein_map <- unique(
  gtf_data[!is.na(gtf_data$gene_id) & !is.na(gtf_data$product), c("gene_id", "product")]
)


pathway_results$Proteins <- sapply(
  strsplit(pathway_results$Genes, ", "), # Split the gene list into individual genes
  function(genes) {
    proteins <- get_proteins_for_genes(genes, gene_protein_map)
    paste(proteins, collapse = ", ") # Combine proteins into a single string
  }
)


# Define a function to process pathways for a drug
process_pathways <- function(go_filtered, go_data, upregulated_genes, gene_protein_map) {
  # Get all enriched pathways (GO IDs)
  enriched_pathways <- go_filtered$GO.ID
  
  # Initialize a data frame to store results
  pathway_results <- data.frame(
    GO_ID = character(),
    Pathway_Name = character(),
    Genes = character(),
    Proteins = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each pathway
  for (pathway_id in enriched_pathways) {
    # Get pathway name/description
    pathway_name <- Term(GOTERM[[pathway_id]])
    
    # Get associated genes
    genes <- intersect(
      genesInTerm(go_data, pathway_id)[[1]],
      upregulated_genes
    )
    
    # Get associated proteins
    proteins <- gene_protein_map$protein_id[gene_protein_map$gene_id %in% genes]
    
    # Combine results into the data frame
    pathway_results <- rbind(
      pathway_results,
      data.frame(
        GO_ID = pathway_id,
        Pathway_Name = pathway_name,
        Genes = paste(genes, collapse = ", "),
        Proteins = paste(unique(proteins), collapse = ", "),
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(pathway_results)
}

# Process two drugs
pathway_results_drug1 <- process_pathways(
  meropenem_up_GO_filtered, meropenem_up_GO_data, 
  meropenem_upregulated_genes, gene_protein_map
)

pathway_results_drug2 <- process_pathways(
  imipenem_up_GO_filtered, imipenem_up_GO_data, 
  imipenem_upregulated_genes, gene_protein_map
)

# Compare results between two drugs
compare_pathways <- function(results1, results2) {
  # Merge on GO_ID to find common pathways
  merged_results <- merge(results1, results2, by = "GO_ID", suffixes = c("_drug1", "_drug2"))
  
  # Find common genes for each pathway
  merged_results$Common_Genes <- mapply(
    function(genes1, genes2) {
      common <- intersect(strsplit(genes1, ", ")[[1]], strsplit(genes2, ", ")[[1]])
      paste(common, collapse = ", ")
    },
    merged_results$Genes_drug1,
    merged_results$Genes_drug2
  )
  
  # Filter out pathways with no common genes
  merged_results <- merged_results[merged_results$Common_Genes != "", ]
  
  return(merged_results)
}

# Perform comparison
common_pathways <- compare_pathways(pathway_results_drug1, pathway_results_drug2)




# Function to process pathways and get shared genes
process_shared_pathways <- function(go_filtered1, go_data1, upregulated_genes1,
                                     go_filtered2, go_data2, upregulated_genes2) {
  # Get all enriched pathways (GO IDs) for both drugs
  pathways1 <- go_filtered1$GO.ID
  pathways2 <- go_filtered2$GO.ID
  
  # Identify common pathways
  common_pathways <- intersect(pathways1, pathways2)
  
  # Initialize a data frame to store results
  shared_results <- data.frame(
    GO_ID = character(),
    Pathway_Name = character(),
    Common_Genes = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each common pathway
  for (pathway_id in common_pathways) {
    # Get pathway name
    pathway_name <- Term(GOTERM[[pathway_id]])
    
    # Get associated genes for both drugs
    genes1 <- intersect(genesInTerm(go_data1, pathway_id)[[1]], upregulated_genes1)
    genes2 <- intersect(genesInTerm(go_data2, pathway_id)[[1]], upregulated_genes2)
    
    # Find common genes
    common_genes <- intersect(genes1, genes2)
    
    # Add to results if common genes exist
    if (length(common_genes) > 0) {
      shared_results <- rbind(
        shared_results,
        data.frame(
          GO_ID = pathway_id,
          Pathway_Name = pathway_name,
          Common_Genes = paste(common_genes, collapse = ", "),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  return(shared_results)
}

# Process and find shared pathways
shared_pathway_results <- process_shared_pathways(
  meropenem_down_GO_filtered, meropenem_down_GO_data, meropenem_downregulated_genes,
  imipenem_down_GO_filtered, imipenem_down_GO_data, 
  imipenem_downregulated_genes
)

# Define the function to map proteins to pathways
map_proteins_to_pathways <- function(pathway_results, gene_protein_map) {
  # Helper function to get unique proteins for a list of genes
  get_proteins_for_genes <- function(genes, map) {
    matched_proteins <- map$product[map$gene_id %in% genes]
    unique(matched_proteins)
  }
  
  # Map proteins for each pathway in the pathway_results data frame
  pathway_results$Proteins <- sapply(
    strsplit(pathway_results$Common_Genes, ", "), # Split the gene list into individual genes
    function(genes) {
      proteins <- get_proteins_for_genes(genes, gene_protein_map)
      paste(proteins, collapse = ", ") # Combine proteins into a single string
    }
  )
  
  return(pathway_results)
}

# Apply the function to shared_pathway_results
final_shared_results <- map_proteins_to_pathways(shared_pathway_results, gene_protein_map)

# Identify overlapping pathways
common_terms <- intersect(imipenem_up_GO_filtered_arranged$Term, meropenem_up_GO_filtered_arranged$Term)

# Filter datasets for overlapping pathways
filtered_data <- imipenem_up_GO_filtered_arranged %>% filter(Term %in% common_terms)


filtered_data %>% 
  ggplot(aes(y = Term, x = GeneRatio)) +  # Swap x and y axes
  geom_col(aes(fill = weight01), width = 0.1) + # Thin bars
  geom_point(aes(size = Significant), color = "black") + # Overlay points
  geom_text(aes(label = Term, x = 0), # Move labels to the left of bars
            hjust = -0.03, vjust = 1.5, size = 5.5, angle = 0) +  # Adjust text angle for horizontal bars
  scale_y_discrete(limits = common_terms) + # Use fixed common_terms
  scale_fill_gradient(low = "red", high = "blue", name = "P-value") + # Gradient fill
  scale_size_continuous(range = c(3, 8), name = "Number of Significant\n         Genes") + # Adjust point size
  theme_light() +
  labs(x = "Enrichment Ratio", # X-axis now shows Enrichment Ratio
       y = "GO Term Description") +  # Y-axis now shows GO Term
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_blank(),  # Adjust y-axis text for better visibility
    axis.text.x = element_text(size = 12),  # Adjust x-axis text for better visibility
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.border = element_rect(color = "black"), 
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    plot.margin = margin(2, 2, 2, 2), # Adjust margins
    aspect.ratio = 1/2 # Maintain aspect ratio
  ) +
  scale_x_continuous(
    limits = c(0, 1),  # Keep the x-axis limit consistent
    breaks = seq(0, 1, 0.25), 
    expand = c(0, 0)
  ) # Consistent x-axis scaling

#Save go gene data
write.table(imipenem_up_genes, file=paste(paste(antibiotic2, "GOupGenes", sep="_"), "tsv", sep="."), quote=FALSE, sep='\t', col.names = NA)
write.table(meropenem_up_genes, file=paste(paste(antibiotic3, "GOupGenes", sep="_"), "tsv", sep="."), quote=FALSE, sep='\t', col.names = NA)
