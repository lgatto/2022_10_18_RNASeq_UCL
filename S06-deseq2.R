#library(rWSBIM2122)
library(tidyverse)
library(DESeq2)

# If files do not exits
#rWSBIM2122::prepare_shell()

load("wsbim2122_data/deseq2/counts.rda")
load("wsbim2122_data/deseq2/coldata.rda")
coldata
head(counts)
dim(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Condition)
dds
assay(dds)[1:5,]
colData(dds)
rowData(dds)

# NB: can also be created from a SE object
se <- SummarizedExperiment(assays = as.matrix(counts),
                           colData = coldata)
dds <- DESeqDataSet(se, design = ~ Condition)


#################################################
# Exercice
#################################################
# Access the count data from the dds object and
# plot the counts distributions of each sample.
as_tibble(assay(dds[, dds$Condition == "KD"])) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>% 
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)


#################################################
# Run DESeq2
# dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

dds <- DESeq(dds)

# rld
rld <- rlogTransformation(dds)
assay(rld)[1:5,]

##################################################
#Pick a gene highly expressed and a gene lowly expressed 
# and inspect the effect of the rlog transformation
##################################################

as_tibble(rowData(dds), rownames = "genes") %>% 
  arrange(baseMean) %>% 
  filter(baseMean > 0) %>% 
  pull(genes) %>% head(1)

as_tibble(assay(dds), rownames = "gene") %>% 
  arrange(desc(sample1))

assay(dds["ENSG00000198804",])  
log2(assay(dds["ENSG00000198804",]))
assay(rld["ENSG00000198804",]) 

assay(dds["ENSG00000248671",])  
log2(assay(dds["ENSG00000248671",]))
assay(rld["ENSG00000248671",]) 
mean(log2(assay(dds["ENSG00000248671",]) + 1))


##################################################


plotPCA(rld, intgroup = "Condition")

pca_data <- plotPCA(rld, intgroup = "Condition", returnData = TRUE)
pca_data

ggplot(pca_data, aes(x = PC1, y = PC2, 
                     shape = Condition, 
                     color = Condition)) +
  geom_point(size = 3) 

ggplot(pca_data, aes(x = PC1, y = PC2, 
                     shape = Condition, 
                     color = Condition,
                     label = name)) +
  geom_point(size = 3) +
  geom_text()

colData(dds)
# Inspecting size factors
sizeFactors(dds)
colData(dds) # stored in the colData






##################################################
# Compare Size Factors to sequencing depth.
##################################################

colSums(assay(dds)) # Compare with sequencing depth


x <- tibble(sample = colnames(dds), 
            SF = sizeFactors(dds),
            SD = colSums(assay(dds)))

SF <- ggplot(x, aes(x = sample, y = SF)) +
  geom_bar(stat = "identity")

SD <- ggplot(x, aes(x = sample, y = SD)) +
  geom_bar(stat = "identity")

library("patchwork")
SF / SD

ggplot(x, aes(x = SD, y = SF)) +
  geom_point() +
  geom_smooth(method = 'lm')




# Dispersion plots
plotDispEsts(dds)


##################################################
# Results
resultsNames(dds)

# Set the mock condition as the ref
dds$Condition <- relevel(dds$Condition, ref = "mock")
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds,
               name = "Condition_KD_vs_mock")

res

res <- results(dds,
        contrast = c("Condition", "KD", "mock"))

results(dds,
        contrast = c("Condition", "mock", "KD"))
               

res_tbl <- as_tibble(res, rownames = "ENSEMBL")
res_tbl
##################################################
# Question 1
# Inspect the results table and identify the 5 “best genes” 
# showing the lowest padjusted value.
##################################################


res_tbl %>%
  arrange(padj) %>%
  head(5) %>% pull(ENSEMBL)

# idem 
head(res_tbl[order(res_tbl$padj), ], 5)

##################################################
# Question 2
# Calculate the mean expression level of these 5 “best genes” 
# using the function counts(). Compare with baseMean values.
##################################################

best_genes <- res_tbl %>%
  arrange(padj) %>%
  head(5) %>% 
  pull(ENSEMBL)

counts(dds[best_genes])
counts(dds[best_genes], normalize = TRUE)

# calculate rowMeans from normalised counts
rowMeans(counts(dds[best_genes], normalize = TRUE))

# => Idem that the Basemean column
res_tbl %>% 
  filter(ENSEMBL %in% best_genes)

##################################################
# Question 3
# Extract the ß coefficient of these 5 “best genes” 
# from the GLM using the function coefficients(). 
# Compare with log2FoldChange values.
##################################################

coefficients(dds)[best_genes, ]
res_tbl %>% 
  filter(ENSEMBL %in% best_genes)

#=> The log2FoldChange corresponds to the Condition_KD_vs_mock coeff

##################################################
# Question 4
# Using the function counts(), evaluate the mean expression 
# levels of these 5 “best genes” in mock cells. 
# Compare with ß coefficients.
##################################################

counts(dds[best_genes, dds$Condition == "mock"])
counts(dds[best_genes, 1:3])
counts(dds[best_genes, dds$Condition == "mock"], normalize = TRUE)
log2(rowMeans(counts(dds[best_genes, dds$Condition == "mock"], normalize = TRUE)))

coefficients(dds[best_genes])









##################################################
# Question 5
# Evaluate the mean expression levels of these 5 “best genes” 
# in KD cells. Compare with ß coefficients.
##################################################
log2(rowMeans(counts(dds[best_genes, dds$Condition == "KD"], normalize = TRUE)))
counts(dds[best_genes, 4:6], normalize = T)

log2(rowMeans(counts(dds[best_genes, 4:6],normalize = T)))
coefficients(dds[best_genes])

log2(mean(counts(dds[best_genes[1], 1:3],normalize = T))) +
  coefficients(dds[best_genes[1]])[,2]

##################################################
# Question 6
# How many genes have no padjusted value? Why?
##################################################

res_tbl %>% 
  filter(is.na(padj)) 



##################################################
# Independant Filtering
metadata(res)$filterThreshold
metadata(res)$filterNumRej

as_tibble(metadata(res)$filterNumRej) %>%
  ggplot(aes(x = theta, y = numRej)) +
  geom_point() +
  geom_vline(xintercept = metadata(res)$filterTheta,
             color = 'red')

##################################################
# Question IF_1
# Actually many of these genes would have been filtered 
# anyway because their basemean == 0. 
# Evaluate how many genes were really filtered by the 
# independent filtering procedure.
##################################################
res_tbl %>% filter(baseMean == 0)


res_tbl %>%
  filter(baseMean > 0 & baseMean < metadata(res)$filterThreshold) 

# idem
res_tbl %>% 
  filter(is.na(padj)) %>% 
  filter(!is.na(pvalue))

##################################################
# Question IF_2
# Re-run the results() function on the same dds object, 
# but set the independent filtering parameter to FALSE. 
# Check how many genes have no padj?
##################################################

res_no_IF <- results(dds, independentFiltering = FALSE)
as_tibble(res_no_IF, rownames = "ENSEMBL") %>% 
  filter(is.na(padj))

# no value for filterThreshold
metadata(res_no_IF)$filterThreshold

##################################################
# Question IF_3
# Imagine another way of filtering genes with very low counts
##################################################

filtering_thr <- 5
# keep genes with counts > 5 in 3 or more samples
keep <- rowSums(counts(dds, normalized = TRUE) >= filtering_thr) >=3
dds_bis <- DESeq(dds[keep, ])

res_bis <- results(dds_bis,
                   name = "Condition_KD_vs_mock",
                   independentFiltering = FALSE)



# No genes without padjusted value
# => logical as we kept only those with at least 5 counts in at least 3 sample
as_tibble(res_bis, rownames = "ENSEMBL") %>%
  filter(is.na(padj)) 
nrow(res_bis)

##################################################

# pvalue histogram
hist(res_tbl$pvalue)

res_tbl %>% 
  filter(pvalue >= 0.8 & pvalue < 0.85)

# the pvalue histogram should be drawn on genes that passed
# the IF procedure
res_tbl %>%
  filter(baseMean > metadata(res)$filterThreshold) %>%
  pull(pvalue) %>%
    hist()

# MA-plot 
plotMA(res)

# Volcano-plot
res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  theme(legend.position = "bottom")

# Count-plot

best_genes <- res_tbl %>%
  arrange(padj)  %>%
  head(4)

as_tibble(counts(dds[best_genes$ENSEMBL, ], normalize = TRUE),
          rownames = 'ENSEMBL') %>%
  pivot_longer(names_to = "sample", values_to = "counts", -ENSEMBL) %>% 
  full_join(as_tibble(coldata, rownames = "sample")) %>% 
  ggplot(aes(x = sample, y = counts, fill = Condition)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free") +
  theme(axis.text.x = element_text(size = 7, angle = 90))


##################################################
# Question
# Identify and inspect counts of the genes plotted in red 
# in the following volcano-plot. These genes have a very 
# large log2FC (|log2FC| > 5) but are far from bearing the 
# lowest padjusted value (their padjusted value is between 0.05 and 1e-5).
##################################################

selected_genes <- res_tbl %>%
  filter(padj < 0.05 & padj > 1e-5 & abs(log2FoldChange) > 5)

selected_genes

as_tibble(counts(dds[selected_genes$ENSEMBL, ], normalize = TRUE),
          rownames = 'ENSEMBL') %>%
  gather(sample, counts, -ENSEMBL) %>%
  left_join(as_tibble(coldata, rownames = "sample")) %>%
  ggplot(aes(x = sample, y = counts, fill = Condition)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ ENSEMBL, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

##################################################
# Question
# Using dispersions() function, compare dispersion 
# values for both group of genes
##################################################

dispersions(dds[selected_genes$ENSEMBL,])
dispersions(dds[best_genes$ENSEMBL,])



##################################################
res_tbl
BiocManager::install(c("biomaRt", "org.Hs.eg.db"))
library("biomaRt")
library("org.Hs.eg.db")

# List of possible datasets can be retrieved 
datasets <- listDatasets(useMart("ensembl"))


## Load homo sapiens ensembl dataset
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mart
## Attributes define the values we are interested to retrieve.
attributes <- listAttributes(mart)
attributes

ensembl_to_geneName <- getBM(attributes = c("ensembl_gene_id", 
                                            "external_gene_name",
                                            "entrezgene_id"),
                             mart = mart)

head(ensembl_to_geneName)
names(ensembl_to_geneName) <- c("ENSEMBL", "gene", "ENTREZID")

res_tbl2 <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj) 

res_tbl2 %>%
  dplyr::select(ENSEMBL, gene, ENTREZID, everything())


  
  
  
  




