library(tidyverse)
source("scripts/R/script_simulation_multivariate-methods.R")

setting = 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure3_4/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
    methods = c("naif", "data.thinning", "Cond"),
    clusterings = c("HAC"),
    nb_experiments = 500,
    setting = setting,
    ns = c(10, 100, 500),
    ps = c(2, 100),
    as = c(0:10),
    Ndraws = c(0, 1000),
    epsilon_dts = c(0.7),
    sigmas = c(1),
    rhos = c(0),
    type_estimation = "known",
    path = path,
    nb_workers = min(40, availableCores()/4)
)

####### Plots #######

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

####### Figure 3 #######
p_H0_ward <- df_clean %>%
  filter(a == 0) %>%
  ggplot(aes(x = pvaleur, linetype = code,
             color = code)) +
  stat_ecdf(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  facet_grid(p ~ n,
             labeller = label_bquote(rows = "m ="~ .(p),
                                     cols = "n ="~ .(n))) +
  labs(x = "p-values", color = "Methods", linetype = "Methods")
# p_H0_ward
ggsave(plot = p_H0_ward, "figure/Figure3.pdf", width = 8, height = 4)


####### Figure 4 #######

df_power <- df_clean %>%
  filter(!(method_inference %in% c("naif"))) %>%
  mutate(reject = pvaleur <= 0.05) %>%
  group_by(n,p,a,ndraws, method_inference, method_clustering, exact, epsilon,
           sigma, rho, method_inference_abre, method_clustering_abre,
           code, code_inf_IS) %>%
  summarise(nb_exp = n(),
            power = mean(reject, na.rm = TRUE),
            mean_time = mean(time, na.rm = TRUE),
            mean_ARI = mean(ARI, na.rm = TRUE))

p_power_ARI <- df_power %>%
  filter(p == 2) %>%
  pivot_longer(cols = c("power", "mean_ARI"), names_to = "categ",
               values_to = "Indicator") %>%
  mutate(categ = recode(categ, "mean_ARI" = "Mean of ARI", .default = categ)) %>%
  mutate(categ = recode(categ, "power" = "Statistical power", .default = categ)) %>%
  mutate(categ = factor(categ, levels = c("Statistical power", "Mean of ARI"))) %>%
  ggplot(aes(x = a, y = Indicator, linetype = code,
             color = code)) +
  geom_point() + geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  labs(x = "a", color = "Methods", y = "",
       linetype = "Methods")+
  facet_grid(categ ~ n,
             labeller = label_bquote(
               cols = "n ="~ .(n)))
# p_power_ARI
ggsave(plot = p_power_ARI, "figures/Figure4.pdf", width = 8, height = 4)
