library(tidyverse)
source("scripts/R/script_simulation_multivariate-methods.R")

setting = 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure5_6/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
  methods = c("data.thinning", "Cond"),
  clusterings = c("kmeans"),
  nb_experiments = 500,
  setting = setting,
  ns = c(10, 100, 200, 500),
  ps = c(2, 50, 100),
  as = c(0:10),
  Ndraws = c(0, 100, 200, 500, 1000),
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


####### Figure 5 #######

df_power <- df_clean %>%
  mutate(reject = pvaleur <= 0.05) %>%
  mutate(ndraws = as.factor(ndraws)) %>%
  group_by(n,p,a,ndraws, method_inference, method_clustering, exact, epsilon,
           sigma, rho, method_inference_abre, method_clustering_abre,
           code, code_inf_IS) %>%
  summarise(nb_exp = n(),
            power = mean(reject, na.rm = TRUE),
            mean_time = mean(time, na.rm = TRUE),
            mean_ARI = mean(ARI, na.rm = TRUE))

p_ARI <- df_power %>%
  filter(p == 2) %>%
  filter(n %in% c(10, 100, 500)) %>%
  pivot_longer(cols = c("power", "mean_ARI"), names_to = "categ",
               values_to = "Indicator") %>%
  mutate(categ = recode(categ, "mean_ARI" = "Mean Adjusted Rand Index", .default = categ)) %>%
  mutate(categ = recode(categ, "power" = "Statistical power", .default = categ)) %>%
  mutate(categ = factor(categ, levels = c("Statistical power", "Mean Adjusted Rand Index"))) %>%
  filter(method_clustering_abre == "KM") %>%
  ggplot(aes(x = a, y = Indicator, linetype = code,
             color = code)) +
  geom_point() + geom_line() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  labs(x = "a", color = "Methods", y = "",
       linetype = "Methods")+
  facet_grid(categ ~ n,
             labeller = label_bquote(
               cols = "n ="~ .(n)))
# p_ARI
ggsave(plot = p_ARI, "figures/Figure5.pdf", width = 8, height = 6)


####### Figure 6 #######

df_power %>%
  filter(method_inference == "Gao") %>%
  filter(a == 0) %>%
  ggplot(aes(x = n, y = mean_time, group = ndraws,
             color = (ndraws))) +
  geom_point() + geom_line() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_y_continuous(transform = "log10") +
  scale_color_manual(values = c("#00BFFF", "#208BB9", "#215B77", "#192F3B", "#000000")) +
  labs(y = "Computation time in sec.",
       color = "Draws")+
  facet_grid( . ~ p,
             labeller = label_bquote(cols = m~"="~ .(p))
  ) -> p_time_kmeans
# p_time_kmeans
ggsave(plot = p_time_kmeans, "figure/Figure6.pdf", width = 3*3, height = 3)

