## BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)
library(airway)
data("airway")
airway

dim(airway)

colData(airway)

rowData(airway)

assay(airway)

assay(airway)[1:5, 1:3]


## get the counts for 
## samples treated with 
## dex and genes 5 419
## and 460

i <- c("ENSG00000000005",
  "ENSG00000000419",
  "ENSG00000000460")

airway[i, ]

airway[1:10, airway$dex == "trt"]

colData(airway)$dex == "trt"

j <- airway$dex == "trt"

airway[i, c(2, 4, 6, 8)]

se2 <- airway[i, j]

assay(se2)

## row ranges

rowRanges(airway)

rowRanges(airway)[[1]]

## build your SE

fls <- dir("wsbim2122_data/count_data/", 
    full.names = TRUE)

tmp <- head(read.delim(fls[1]), 10)

View(tmp)

dim(read.delim(fls[1]))

m <- sapply(fls, 
                 function(f) 
                   read.delim(f)[, 7])

class(counts)
dim(counts)
m[1:5, 1:3]

colnames(m) <- paste0("sample", 1:6)
rownames(m) <- read.delim(fls[1])[, 1]

cd <- data.frame(sample = 1:6, 
                 group = rep(c("A", "B"),
                             each = 3),
                 row.names = colnames(counts))
cd

colnames(counts)

## cd <- read.csv("myColData.csv")

se <- SummarizedExperiment(
  assays = list(counts = m),
  colData = cd)

se
se <- SummarizedExperiment(
  assays = list(counts = m),
  colData = cd, 
  rowData = rd)

colData(se)
se

## ~ group

rowData(se)

