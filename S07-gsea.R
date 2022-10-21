library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

## data/results from previous chapter
dds <- readRDS("data/dds.rds")
ensembl_to_geneNames <- readRDS("data/ensembl_to_geneName.rds")
res_tbl <- readRDS("data/res_tbl.rds")

View(res_tbl)

table(res_tbl$padj < 0.05)

library(GO.db)
library(org.Hs.eg.db)

GO_0005925 <- AnnotationDbi::select(
  org.Hs.eg.db,
  key = "GO:0005925",
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "GO") %>% 
  as_tibble() %>% 
  filter(!duplicated(ENTREZID))

View(GO_0005925)

## ORA

de_genes <- res_tbl %>% 
  filter(padj < 0.01,
         abs(log2FoldChange) > 1) %>% 
  pull(ENTREZID)

not_de_genes <- res_tbl[res_tbl$padj > 0.01 | 
  abs(res_tbl$log2FoldChange) < 1, 
  "ENTREZID"] %>% 
  pull(1)
         

length(intersect(res_tbl$ENTREZID[res_tbl$padj < 0.05],
                 GO_0005925$ENTREZID))

n <- length(intersect(de_genes, GO_0005925$ENTREZID))
m <- length(intersect(not_de_genes, GO_0005925$ENTREZID))
p <- length(setdiff(de_genes, GO_0005925$ENTREZID))
q <- length(setdiff(not_de_genes, GO_0005925$ENTREZID))


cont_mat <- matrix(c(n, m, p, q), nrow = 2)
rownames(cont_mat) <- c("DE", "not_DE")
colnames(cont_mat) <- c("GO", "not_GO")
cont_mat

fisher.test(cont_mat, alternative = "greater")
fisher.test(cont_mat, alternative = "two.sided")

## GSEA

res_tbl$inGO <- res_tbl$ENTREZID %in% GO_0005925$ENTREZID

dplyr::select(res_tbl, 
              ENTREZID, padj, inGO) %>% 
  View()

table(res_tbl$inGO)


colData(dds)
sample(colData(dds)$Condition)
combn(6, 3)

resultsNames(dds)
results(dds, name = "Condition_KD_vs_mock") %>% 
  as_tibble() %>% 
  arrange(padj)
  

head(res_tbl)


library(clusterProfiler)

## ORA

go_ora <- enrichGO(gene = de_genes,
         universe = as.character(res_tbl$ENTREZID),
         OrgDb = org.Hs.eg.db,
         ont = "CC",
         pvalueCutoff = 0.1,
         readable = TRUE) %>% 
  as_tibble()

View(go_ora)

ordered_genes <- res_tbl %>% 
  filter(!is.na(ENTREZID)) %>% 
  arrange(desc(abs(stat))) %>% 
  pull(ENTREZID)

ordered_genes <- abs(res_tbl$stat)
names(ordered_genes) <- res_tbl$ENTREZID
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Hs.eg.db,
                 scoreType = "pos") %>%
  as_tibble


View(go_gsea)

## enrichKEGG
## gseKEGG

