suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))

library(stringr) 

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

# Use enums in your code
antibiotic1 <- antibiotics$CIPROFLOXACIN
antibiotic2 <- antibiotics$MOXIFLOXACIN
antibiotic3 <- antibiotics$TOBRAMYCIN

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