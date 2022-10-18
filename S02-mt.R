library(rWSBIM1322)
data(tdata1)
head(tdata1)
dim(tdata1)
class(tdata1)

log_tdata1 <- log2(tdata1)

head(log_tdata1)

limma::plotDensities(tdata1, 
                     legend = FALSE)
limma::plotDensities(log_tdata1, 
                     legend = FALSE)

x <- log_tdata1[73, ]
x
t.test(x[1:3], x[4:6])

t.test(x[1:3], x[4:6])$p.value

my_t_test <- function(xin) {
  t.test(xin[1:3], xin[4:6])$p.value
}

my_t_test(x)

my_t_test(log_tdata1[12, ])

pv <- apply(log_tdata1, 1, my_t_test)
pv ## 100 p-values, one for each feature

## any significant (at alpha 0.05)?

table(pv <= 0.05)

which(pv <= 0.05)

pv[pv <= 0.05]

log_tdata1[which(pv <= 0.05)
, ]

hist(rnorm(10000, 10))

t.test(rnorm(3, 10), rnorm(3, 10))$p.value

adjpv <- p.adjust(pv, method = "BH")

table(adjpv <= 0.05)

hist(pv)
head(pv)
head(adjpv)


hist(pv, breaks = 20)
