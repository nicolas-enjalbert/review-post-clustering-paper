library(tidyverse)
library(latex2exp)
source("scripts/R/script_simulation_multivariate-methods.R")

setting = 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure7_17/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
  methods = c("data.thinning"),
  clusterings = c("HAC"),
  nb_experiments = 500,
  setting = setting,
  ns = c(100),
  ps = c(2, 50, 100),
  as = c(0:10),
  Ndraws = c(0),
  epsilon_dts = c(0.1, 0.3, 0.5, 0.7, 0.9),
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
df_clean = data_cleaned$df_clean %>%
  mutate(epsilon = as.factor(epsilon))
color_palette = data_cleaned$color_palette
linetype_palette = data_cleaned$linetype_palette
color_palette_inf_IS = data_cleaned$color_palette_inf_IS
linetype_palette_inf_IS = data_cleaned$linetype_palette_inf_IS

####### Figure 7 #######

df_clean %>%
  mutate(reject = pvaleur <= 0.05) %>%
  group_by(n,p,a,ndraws, method_inference, method_clustering, exact, epsilon) %>%
  summarise(nb_exp = n(),
            power = mean(reject, na.rm = TRUE),
            mean_time = mean(time, na.rm = TRUE),
            mean_ARI = mean(ARI, na.rm = TRUE)) %>%
  pivot_longer(cols = c("power", "mean_ARI"), names_to = "categ",
               values_to = "Indicator") %>%
  mutate(categ = recode(categ, "mean_ARI" = "Mean Adjusted Rand Index",
                        .default = categ)) %>%
  mutate(categ = recode(categ, "power" = "Statistical power",
                        .default = categ)) %>%
  mutate(categ = factor(categ, levels = c("Statistical power",
                                          "Mean Adjusted Rand Index"))) %>%
  ggplot(aes(x = a, y = Indicator,  color = epsilon)) +
  geom_point() + geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(y = "", color = TeX("$\\epsilon$")) +
  facet_grid(categ ~ p,
             labeller = label_bquote(cols = "m ="~ .(p))) -> p_ARI_dt
ggsave(plot = p_ARI_dt, "figures/Figure7.pdf", width = 8, height = 6)

####### Figure 17 #######

df_clean %>%
  filter(a == 0) %>%
  ggplot(aes(x = pvaleur, color = epsilon)) +
  stat_ecdf() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  facet_grid(n ~ p,
             labeller = label_bquote(cols = m~"="~ .(p),
                                     rows = n~"="~ .(n))) +
  labs(x = "p-values", color = TeX("$\\epsilon$")) -> p_dt_H0
ggsave(plot = p_dt_H0, "figure/Figure17.pdf", width = 9, height = 4)
