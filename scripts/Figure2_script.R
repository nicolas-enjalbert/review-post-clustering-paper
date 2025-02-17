library(tidyverse)
library(ggpubr)
source("scripts/R/utils.R")

############ Setting 1 ##############
n = 500
p = 2
a = 5

{
  signal <- signal_2clusters(n = n, p = p, a = a, prob_clusters = c(0.5, 0.5),
                             nb_variable_signal = p/2, verbose = FALSE)

  MatMean <- signal$X
  truthS <- signal$truthS
  truthF <- signal$truthF

  K <- length(unique(truthS))


  Sigma <- diag(p)


}

X2 <- SimuMultiCorr(n = n, MatMean = MatMean, Sigma = Sigma)

data.frame(X2, Truth = as.factor(truthS)) %>%
  ggplot(aes(x = X1, y = X2, color = Truth)) +
  geom_point(size = 1.9) +
  labs(x = "Dim 1", y = "Dim 2", color = TeX("Labels of the true partition")) +
  theme_bw()  +
  theme(legend.position = "bottom") -> p_setting1
ggExtra::ggMarginal(p_setting1, type = "density", groupColour = TRUE,
                    groupFill = TRUE) -> p_setting1
# p_setting1
ggsave(plot = p_setting1, "figures/Figure2.pdf", width = 9, height = 5)
