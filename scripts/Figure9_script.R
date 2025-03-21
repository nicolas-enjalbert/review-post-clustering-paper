library(tidyverse)
library(ggpubr)
source("scripts/R/utils.R")
set.seed(13012024)
############ Setting 2 ##############
n <- 500
p <- 2
a <- 5

signal <- signal_3clusters_eq(n,
  a = a, p = p,
  prob_clusters = c(1 / 3, 1 / 3, 1 / 3),
  verbose = TRUE
)

MatMean <- signal$X
truthS <- signal$truthS
truthF <- signal$truthF

K <- length(unique(truthS))


Sigma <- diag(p)

X3 <- SimuMultiCorr(n = n, MatMean = MatMean, Sigma = Sigma)

p_setting2 <- data.frame(X3, Truth = as.factor(truthS)) %>%
  ggplot(aes(x = X1, y = X2, color = Truth)) +
  geom_point(size = 1.9) +
  labs(x = "Variable j = 1", y = "Variable j = 2", color = "True partition") +
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "left")
p_setting2 <- ggExtra::ggMarginal(p_setting2,
  type = "density", groupColour = TRUE,
  groupFill = TRUE
)
ggsave(plot = p_setting2, "figures/Figure9.pdf", width = 8, height = 4)
