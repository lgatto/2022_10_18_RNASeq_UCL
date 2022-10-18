library(tidyverse)

counts_result <- read_tsv(file = "wsbim2122_data/count_data/sample1_counts.tsv.gz")

## number of genes

nrow(counts_result)

## How many reads were assigned to genes?

sum(counts_result[[7]])

counts_result <- counts_result %>%
  dplyr::rename(counts = "../processed_data/bam/sample1.bam")
sum(counts_result$counts)

## Inspect the counts. What is the maximum value? 

summary(counts_result[[7]])

range(counts_result[[7]])

min(counts_result[[7]])

max(counts_result[[7]])

## How many genes have zero counts ?

table(counts_result[[7]] == 0)

sum(counts_result[[7]] == 0)

counts_result %>%
  filter(counts == 0) %>% 
  nrow()
  
## Plot the distribution of counts. How does it look like?

counts_result %>%
  ggplot(aes(x = log2(counts + 1))) +
  geom_histogram(bins = 20)

hist(log2(counts_result[[7]] + 1))

counts_result %>%
  ggplot(aes(x = log1p(counts))) +
  geom_histogram(bins = 20)

hist(log1p(counts_result[[7]]))

