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

