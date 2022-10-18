# Course material and installation instructions:
#   
#   https://uclouvain-cbio.github.io/
#   bioinfo-training-02-rnaseq/
#   
#   Start with
# 
# install.packages("remotes")
# 
# and then instructions.

install.packages("remotes")
BiocManager::install("UCLouvain-CBIO/rWSBIM2122")
BiocManager::install("UCLouvain-CBIO/rWSBIM1322")

rWSBIM2122::prepare_shell()

library(rWSBIM2122)
data(gexp)
gexp

t.test(expression ~ group, 
       data = gexp)

library(ggplot2)

p <- ggplot(gexp, 
       aes(x = group, y = expression))  +
  geom_boxplot() + geom_jitter()
p

mod1 <- lm(expression ~ group, 
   data = gexp)
mod1

4.33 - 2.75

p  +
  geom_point(data = data.frame(x = c(1, 2), y = c(2.75, 4.33)),
             aes(x = x, y = y),
             colour = "red", size = 5) +
  geom_abline(intercept = coefficients(mod1)[1] - coefficients(mod1)[2],
              slope = coefficients(mod1)[2])


t.test(expression ~ group, 
       data = gexp, 
       paired = TRUE)

gexp2 <- dplyr::full_join(gexp[1:10, ], gexp[11:20, ],
                          by = "ID", suffix = c("_1", "_2"))

ggplot(gexp, aes(x = group, y = expression)) +
  geom_boxplot() +
  geom_point(aes(colour = ID), size = 4) +
  geom_segment(data = gexp2,
               aes(x = group_1, xend = group_2,
                   y = expression_1, yend = expression_2,
                   colour = ID))
gexp

mod2 <- lm(expression ~ group + ID, 
           data = gexp)
mod2

summary(mod2)