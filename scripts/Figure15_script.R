library(tidyverse)
source("scripts/R/script_simulation_multivariate-methods.R")

setting = 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure15/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
  methods = c("data.thinning", "naif"),
  clusterings = c("HAC", "kmeans", "GMM"),
  nb_experiments = 500,
  setting = setting,
  ns = c(10, 100, 500),
  ps = c(2, 10, 50),
  as = c(0:10),
  Ndraws = c(0),
  epsilon_dts = c(0.7),
  sigmas = c(1),
  rhos = c(0),
  type_estimation = "known",
  path = path,
  nb_workers = min(40, availableCores()/4)
)

####### Figures #######

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

####### Figure 15 #######

df_power <- df_clean %>%
  mutate(true_clustering = ARI == 1) %>%
  group_by(n,p,a,ndraws, method_inference, method_clustering, exact, epsilon,
           sigma, struc_cov, method_inference_abre, method_clustering_abre,
           code, code_inf_IS) %>%
  summarise(nb_exp = n(),
            sum_true_clustering = sum(true_clustering),
            prop_true_clustering = mean(true_clustering))

p_prop_recover_clustering <- df_power %>%
  mutate(code_inf_IS = recode(code_inf_IS,
                              "Naive" = "Entire dataset",
                              "DT" = "Thinned dataset")) %>%
  ggplot(aes(x = a, y = prop_true_clustering,
             linetype = method_clustering_abre,
             color = code_inf_IS)) +
  geom_point() + geom_line() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  facet_grid(n ~ p,
             labeller = label_bquote(rows = "n ="~ .(n),
                                     cols = "m ="~ .(p))) +
  labs(x = "a", y = "Proportion of recovered clustering [through 500 experiments]",
       color = 'Dataset used', linetype = "Clustering method")
ggsave(filename = "figure/Figure15.pdf",
       plot = p_prop_recover_clustering,
       height = 6, width = 8)
