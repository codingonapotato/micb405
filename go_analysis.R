suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))

#### Load GO mapping file and create character vector with all gene IDs from the file ####
geneID2GO <- readMappings("go_mapping.tsv")
geneUniverse <- names(geneID2GO)

#### Load differential expression analysis data ####
colistin_vs_control_90 <- read_csv("de_seq_results/colistin_vs_ctrl_90.csv")
imipenem_vs_control_90 <- read_csv("de_seq_results/imipenem_vs_ctrl_90.csv")
meropenem_vs_control_90 <- read_csv("de_seq_results/meropenem_vs_ctrl_90.csv")

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

#### Build GOdata objects for up/down-regulated data ####
colistin_up_GO_data <- new("topGOdata", 
                  description = "Colisin_90min", 
                  ontology = "BP", 
                  allGenes = colistin_up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)
imipenem_up_GO_data <- new("topGOdata", 
                           description = "Imipenem_90min", 
                           ontology = "BP", 
                           allGenes = imipenem_up_gene_list,
                           annot = annFUN.gene2GO,
                           gene2GO = geneID2GO)
meropenem_up_GO_data <- new("topGOdata", 
                           description = "Meropenem_90min", 
                           ontology = "BP", 
                           allGenes = meropenem_up_gene_list,
                           annot = annFUN.gene2GO,
                           gene2GO = geneID2GO)

colistin_down_GO_data <- new("topGOdata",
                    description = "Colisin_90min",
                    ontology = "BP",
                    allGenes = colistin_down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
imipenem_down_GO_data <- new("topGOdata",
                             description = "Imipenem_90min",
                             ontology = "BP",
                             allGenes = imipenem_down_gene_list,
                             annot = annFUN.gene2GO,
                             gene2GO = geneID2GO)
meropenem_down_GO_data <- new("topGOdata",
                             description = "Meropenem_90min",
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

#### Extract a summary of up/down-regulated gene results ####
colistin_up_GO <- GenTable(colistin_up_GO_data,
                           weight01 = colistin_up_result,
                           orderBy = "colistin_up_result",
                           ranksOf = "colistin_up_result",
                           topNodes = 50)
imipenem_up_GO <- GenTable(imipenem_up_GO_data,
                           weight01 = imipenem_up_result,
                           orderBy = "imipenem_up_result",
                           ranksOf = "imipenem_up_result",
                           topNodes = 50)
meropenem_up_GO <- GenTable(meropenem_up_GO_data,
                           weight01 = meropenem_up_result,
                           orderBy = "meropenem_up_result",
                           ranksOf = "meropenem_up_result",
                           topNodes = 50)

colistin_down_GO <- GenTable(colistin_down_GO_data,
                             weight01 = colistin_down_result,
                             orderBy = "colistin_down_result",
                             ranksOf = "colistin_down_result",
                             topNodes = 50)
imipenem_down_GO <- GenTable(imipenem_down_GO_data,
                             weight01 = imipenem_down_result,
                             orderBy = "imipenem_down_result",
                             ranksOf = "imipenem_down_result",
                             topNodes = 50)
meropenem_down_GO <- GenTable(meropenem_down_GO_data,
                             weight01 = meropenem_down_result,
                             orderBy = "meropenem_down_result",
                             ranksOf = "meropenem_down_result",
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
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = colistin_up_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")
gg_imipenem_up <- imipenem_up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = imipenem_up_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")
gg_meropenem_up <- meropenem_up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = meropenem_up_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")

gg_colistin_down <- colistin_down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = colistin_down_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")
gg_imipenem_down <- imipenem_down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = imipenem_down_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")
gg_meropenem_down <- meropenem_down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = meropenem_down_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")

#### Save visualizations to file ####
ggsave(filename="figures/colistin_upregulated.png", plot = gg_colistin_up,
       width = 7, height = 7)
ggsave(filename="figures/imipenem_upregulated.png", plot = gg_imipenem_up,
       width = 7, height = 7)
ggsave(filename="figures/meropenem_upregulated.png", plot = gg_meropenem_up,
       width = 7, height = 7)

ggsave(filename="figures/colistin_downregulated.png", plot = gg_colistin_down,
       width = 7, height = 7)
ggsave(filename="figures/imipenem_downregulated.png", plot = gg_imipenem_down,
       width = 7, height = 7)
ggsave(filename="figures/meropenem_downregulated.png", plot = gg_meropenem_down,
       width = 7, height = 7)