library(tidyverse)
source("scripts/R/script_simulation_multivariate-methods.R")

setting = 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure8_9/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
  methods = c("data.thinning", "Cond", "YFB"),
  clusterings = c("HAC"),
  nb_experiments = 500,
  setting = setting,
  ns = c(500),
  ps = c(2, 10),
  as = c(0:10),
  Ndraws = c(0, 1000),
  epsilon_dts = c(0.7),
  sigmas = c(1),
  rhos = c(0, 0.3, 0.5),
  type_estimation = c("known", "all", "intra"),
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

data_cleaned <- clean_data_estimation(df_end)
df_clean = data_cleaned$df_clean
color_palette = data_cleaned$color_palette
linetype_palette = data_cleaned$linetype_palette

####### Figure 8 #######

p_H0_dependence <- df_clean %>%
  filter(a == 0) %>%
  ggplot(aes(x = pvaleur, color = code, linetype = code)) +
  stat_ecdf(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(p ~ rho,
             labeller = label_bquote(cols = rho~"="~ .(rho),
                                     rows = m~"="~ .(p))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  labs(x = "p-values", color = "Methods", linetype = "Methods")
p_H0_dependence
ggsave(plot = p_H0_dependence,
       "figures/Figure8.pdf", width = 8, height = 4)

####### Figure 9 #######

df_clean %>%
  filter(!(code == "Cond-intra" & (rho != 0 | p != 10))) %>%
  filter(!(code == "YFB" & rho != 0)) %>%
  filter(!(code == "YFB-IS" & (rho != 0 | p != 2))) %>%
  filter(code != "Cond-intra") %>%
  mutate(reject = pvaleur <= 0.05) %>%
  group_by(n,p,a,ndraws, method_inference, method_clustering, exact, sigma,
           rho, type_estimation, code ) %>%
  summarise(nb_exp = n(),
            power = mean(reject, na.rm = TRUE),
            mean_time = mean(time, na.rm = TRUE),
            mean_ARI = mean(ARI, na.rm = TRUE)) %>%
  ggplot(aes(x = a, y = power, color = code, linetype = code)) +
  geom_point() + geom_line() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette)+
  labs(color = "Methods", linetype = "Methods", y = "Statistical power")+
  facet_grid(p ~ rho,
             labeller = label_bquote(cols = rho~"="~ .(rho),
                                     rows = m~"="~ .(p))
  ) -> p_power_dependence
# p_power_dependence
ggsave(plot = p_power_dependence,
       "figures/Figure9.pdf",
       width = 8, height = 4)
