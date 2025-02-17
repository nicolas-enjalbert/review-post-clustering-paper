library(tidyverse)
library(clusterpval)
library(KmeansInference)
library(mclust)
library(Rmixmod)
library(Hotelling)
source("scripts/R/function_yun.R")
source("scripts/R/utils.R")
library(future.apply)


path <- sprintf("results/simulation_multivariate_setting%s/Figure18_19/", setting)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

#####################

nb_workers <- min(40, availableCores()/2)

plan(multisession, workers = nb_workers)
set.seed(130124)

nb_experiments = 500

# Possible values for the number of sample
ns <- c(30, 500)
# Possible values for the number of variable
ps <- c(2, 100)
# Possible values for the signal a = || \eta^T X ||_2 (independent from p)
as <- c(0)

# value of sigma
sigmas <- c(1)
rhos = 0
configs <- expand.grid(ns = ns,
                       ps = ps,
                       as = as,
                       sigmas = sigmas,
                       rhos = rhos)
seq_configs <- 1:nrow(configs)

K = 2
HACMethods <- c("ward.D")
ndraws <- 1000 #Monte carlo sample
alphas <- c(0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 1)
cc = 1
for(cc in seq_configs){
  config <- configs[cc,]
  a <- config[["as"]]
  n <- config[["ns"]]
  p <- config[["ps"]]
  sigma <- config[["sigmas"]]
  rho <- config[["rhos"]]
  iso = rho == 0

  print(sprintf("a = %s ; n = %s ; p = %s ; sigma = %s ; rho = %s",
                a, n, p, sigma, rho))


  res <- future.apply::future_lapply(1:nb_experiments, future.seed = TRUE, FUN = function(nexp) {

    {
      signal <- signal_2clusters(n = n, p = p, a = a,
                                 prob_clusters = c(0.5, 0.5),
                                 nb_variable_signal = p/2, verbose = FALSE)

      MatMean <- signal$X
      truthS <- signal$truthS
      truthF <- signal$truthF

      K <- length(unique(truthS))

      sigmas <- rep(1, p)
      Sigma <- diag(p)
      SigmaInv <- diag(p) # solve(Sigma) in the general case


    }

    X <- SimuMultiCorr(n = n, MatMean = MatMean, Sigma = Sigma)

    df <- data.frame()

    ####### yun HAC ward

    for(alpha in alphas){
      t1 <- Sys.time()
      yun_ISMC = fun_proposed_approx(X, K = K, k1 = 1, k2 = 2, ndraws = ndraws,
                                     alpha = alpha, cl_meth = "ward")
      time_ <- difftime(Sys.time(), t1, units = "secs")
      ARI <- mclust::adjustedRandIndex(yun_ISMC$cl, truthS)
      df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = ndraws,
                                 a = a, acc_rate = yun_ISMC$acc_rate,
                                 alpha = alpha,
                                 method_inference = "unknown_variance",
                                 method_clustering = "HAC_ward.D",
                                 exact = FALSE,
                                 sigma = sigma,  rho = rho,
                                 pvaleur = yun_ISMC$pval,
                                 stat = 0,
                                 ARI = ARI, time = as.numeric(time_ )))

    }



    return(df)
  }
  )
  df_end <- Reduce(rbind, (res))

  name_file <- sprintf("sim_multi_set1_K=%s_a=%s_n=%s_p=%s_nb_experiment=%s.RDS",
                       K, a, n, p, nb_experiments)
  file.path <- paste(path, name_file, sep = "/")
  saveRDS(df_end, file = file.path)
}


######################################################
K=2
pattern <- sprintf("sim_multi_set1_K=%s_a=(.*)_n=(.*)_p=(.*)_nb_experiment=(.*).RDS",
                   K)
filenames <- list.files(path, pattern = pattern)

df_end <- data.frame()
for(filename in filenames){
  df <- readRDS(paste(path, filename, sep = "/"))
  df_end <- rbind(df_end, df)
}



p_H0 <- df_H0 %>%
  filter(a == 0) %>%
  mutate(alpha = as.factor(alpha))%>%
  ggplot(aes(x = pvaleur, color = alpha)) +
  stat_ecdf(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  facet_grid(p ~ n,
             labeller = label_bquote(rows = m~"="~ .(p),
                                     cols = n~"="~ .(n))) +
  labs(x = "p-values", linetype = "Computation of the p-value", color = "gamma")
p_H0
ggsave(plot = p_H0, "figures/Figure18.pdf", height = 5, width = 7)




df_end %>%
  filter(a == 0) %>%
  mutate(alpha = as.factor(alpha)) %>%
  ggplot(aes(y = acc_rate, x = alpha, color = alpha, group = alpha)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_grid(p ~ n,
             labeller = label_bquote(rows = m~"="~ .(p),
                                     cols = n~"="~ .(n))) +
  labs(y = "Accuracy rate", x = "gamma", color = "gamma") -> p_boxplot_alpha
p_boxplot_alpha
ggsave(plot = p_boxplot_alpha,
       filename = "figures/Figure19.pdf",
       height = 5,
       width = 7)
