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
