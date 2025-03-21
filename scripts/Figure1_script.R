library(tidyverse)
library(ggpubr)
library(latex2exp)
source("scripts/R/utils.R")

set.seed(13012026)
n <- 500
p <- 2
a <- 0

signal <- signal_2clusters(
  n = n, p = p, a = a, prob_clusters = c(0.5, 0.5),
  nb_variable_signal = p / 2, verbose = FALSE
)

MatMean <- signal$X
truthS <- signal$truthS
truthF <- signal$truthF

K <- length(unique(truthS))

Sigma <- diag(rep(1, p))

X2 <- SimuMultiCorr(n = n, MatMean = MatMean, Sigma = Sigma)
cl <- cutree(hclust(dist(X2, method = "euclidean"), method = "ward.D2"), 2)

## plot 1
p_setting1 <- data.frame(X2, Truth = as.factor(cl)) %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(size = 1.9, color = "darkgrey") +
  labs(x = "Variable j = 1", y = "Variable j = 2") +
  geom_point(x = 0, y = 0, color = "black") +
  theme(legend.position = "bottom")
p_setting1 <- ggExtra::ggMarginal(p_setting1, type = "density", fill = "grey")
# p_setting1

## plot 2
cl1 <- colMeans(X2[cl == 1, ])
cl2 <- colMeans(X2[cl == 2, ])
p_setting1_clustering <- data.frame(X2, Truth = as.factor(cl)) %>%
  ggplot(aes(x = X1, y = X2, color = Truth)) +
  geom_point(size = 1.9) +
  labs(x = "Variable j = 1", y = "Variable j = 2", color = "Clusters (HAC)") +
  geom_point(x = cl1[1], y = cl1[2], color = "black") +
  geom_point(x = cl2[1], y = cl2[2], color = "black") +
  theme(
    legend.position = c(0.85, 0.15),
    legend.text = element_text(size = 6), # Taille du texte
    legend.title = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p_setting1_clustering <- ggExtra::ggMarginal(p_setting1_clustering,
  type = "density", groupColour = TRUE,
  groupFill = TRUE
)
# p_setting1_clustering


## plot3
## simulation of 500 dataset and use the Wald test
nb_experiments <- 500
pvalues <- c()
for (i in seq(nb_experiments)) {
  signal <- signal_2clusters(
    n = n, p = p, a = a, prob_clusters = c(0.5, 0.5),
    nb_variable_signal = p / 2, verbose = FALSE
  )

  MatMean <- signal$X
  truthS <- signal$truthS
  truthF <- signal$truthF

  K <- length(unique(truthS))
  sigma <- 1
  Sigma <- diag(rep(sigma, p))


  X2 <- SimuMultiCorr(n = n, MatMean = MatMean, Sigma = Sigma)
  cl <- cutree(hclust(dist(X2, method = "euclidean"), method = "ward.D2"),
    k = K
  )

  eta <- eta_comparison_mean(cl, 1, 2)
  pvalues[i] <- multivariate_Z_test(X2, eta, sigma)$pval
}

p_H0_ward <- data.frame(pvalues = pvalues, code = "Naive-HAC") %>%
  ggplot(aes(x = pvalues, color = code)) +
  stat_ecdf(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0, 1) +
  theme(
    legend.position = c(0.85, 0.1),
    legend.box = "vertical",
    legend.text = element_text(size = 6), # Taille du texte
    legend.title = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = "#C77CFF") +
  labs(
    x = "p-values", color = "Methods", linetype = "Methods",
    y = "ECDF (through 500 experiments)"
  )

p_naive <- ggarrange(p_setting1, p_setting1_clustering, p_H0_ward,
  ncol = 3, nrow = 1
)
ggsave(
  plot = p_naive, "figures/Figure1.pdf",
  width = 12, height = 4
)
