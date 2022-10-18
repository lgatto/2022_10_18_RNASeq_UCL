library(pRolocdata)
data(mulvey2015_se)

mulvey2015_se

colData(mulvey2015_se)
dim(mulvey2015_se)


x <- assay(mulvey2015_se)
dim(x)
pca <- prcomp(t(x), scale = TRUE)

pca
class(pca)
str(pca)

head(pca$x)

plot(pca$x[, 1:2])
plot(pca$x[, c(1, 3)])
