library(tidyverse)
source("scripts/R/script_simulation_multivariate-methods.R")

setting = 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure14/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
  methods = c("naif", "data.thinning", "Cond"),
  clusterings = c("HAC", "kmeans", "GMM"),
  nb_experiments = 500,
  setting = setting,
  ns = c(10, 100, 500),
  ps = c(2),
  as = c(0),
  Ndraws = c(0, 1000),
  epsilon_dts = c(0.7),
  sigmas = c(1),
  rhos = c(0),
  type_estimation = "known",
  path = path,
  nb_workers = min(40, availableCores()/4)
)

####### Figure #######

pattern <- "sim_multi_set1_K=2_a=(.*)_n=(.*)_p=(.*)_nb_experiment=(.*)_sigma=(.*)_rho=(.*).RDS"
filenames <- list.files(path, pattern = pattern)

df_end <- data.frame()
for(filename in filenames){
  df <- readRDS(paste(path, filename, sep = "/"))
  df_end <- rbind(df_end, df)
}

data_cleaned <- clean_data(df_end)
df_clean = data_cleaned$df_clean
color_palette = data_cleaned$color_palette
linetype_palette = data_cleaned$linetype_palette
color_palette_inf_IS = data_cleaned$color_palette_inf_IS
linetype_palette_inf_IS = data_cleaned$linetype_palette_inf_IS

####### Figure 14 #######
p_H0 <- df_clean %>%
  ggplot(aes(x = pvaleur, linetype = code_inf_IS,
             color = code_inf_IS)) +
  stat_ecdf(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical")  +
  facet_grid(method_clustering ~ n,
             labeller = label_bquote(
               cols = n~"="~ .(n))) +
  labs(x = "p-values", color = "Methods", linetype = "Methods") +
  scale_color_manual(values = color_palette_inf_IS) +
  scale_linetype_manual(values = linetype_palette_inf_IS)
# p_H0
ggsave(plot = p_H0, "figures/Figure14.pdf", width = 8, height = 6)

