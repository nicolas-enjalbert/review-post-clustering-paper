library(tidyverse)
source("scripts/R/script_simulation_marginal_methods.R")

setting = 2
path <- sprintf("results/simulation_marginal_setting%s/Figure11_12_13/", setting)

####### Simulations #######
comparison_univariate_inference_methods(
    methods = c("Z.test", "data.thinning", "Hivert", "Chen", "tn.test", "poclin"),
    clusterings = c("HAC"),
    nb_experiments = 500,
    setting = setting,
    ns = 100,
    ps = 3,
    as = c(0:10),
    Ndraws = c(0, 1000),
    epsilon_dts = c(0.7),
    sigmas = c(1),
    rhos = c(0),
    path = path,
    nb_workers = min(40, availableCores()/4)
)

####### Figures #######

pattern <- "sim_univariate_set2_K=3_a=(.*)_n=(.*)_p=(.*)_nb_experiment=(.*)_sigma=(.*)_rho=(.*).RDS"
filenames <- list.files(path, pattern = pattern)

df_end <- data.frame()
for(filename in filenames){
  df <- readRDS(paste(path, filename, sep = "/"))
  df_end <- rbind(df_end, df)
}

data_cleaned <- clean_data_marginal(df_end)
df_clean = data_cleaned$df_clean
color_palette = data_cleaned$color_palette

####### Figure 11 #######

p_H0 <- df_clean %>%
  filter(a %in% c(0,3,5)) %>%
  filter(method_inference != "Hivert - direct_test - estimated") %>%
  filter( H0 %in% c("H0 j = 2", "H0 j = 3")) %>%
  mutate(H0 = recode(H0,
                     "H0 j = 2" = "j = 2",
                     "H0 j = 3" = "j = 3")) %>%
  ggplot(aes(x = pvalues, color = code)) +
  stat_ecdf(linewidth = 1.2) +
  facet_grid(H0 ~ a,
             labeller = label_bquote(
               rows = .(H0) * "   (under " * H[0] * ")",
               cols = a: .(a))
  ) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "p-values", color = "Methods")+
  scale_color_manual(values = color_palette)
# p_H0
ggsave(plot = p_H0,
       filename = "figures/Figure11.pdf",
       width = 8, height = 6)

####### Figure 12 #######


df_power <- df_clean %>%
  mutate(reject = pvalues <= 0.05) %>%
  group_by(n,p,a,ndraws, method_inference, method_clustering, exact, epsilon,
           sigma, rho, H0, variables, truesignal, code, TV) %>%
  summarise(nb_exp = n(),
            power = mean(reject, na.rm = TRUE),
            mean_time = mean(time, na.rm = TRUE),
            mean_ARI = mean(ARI, na.rm = TRUE))


p_ARI <- df_power %>%
  filter(variables == 3) %>%
  filter(n == 100) %>%
  filter(method_inference != "Hivert - direct_test - estimated") %>%
  ggplot(aes(x = a, y = mean_ARI, color = code)) +
  geom_point() + geom_line() +
  # facet_wrap(vars(n)) +
  scale_color_manual(values = color_palette)  +
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(y = "Average of Adjusted Rand Index (through 500 experiments)", color = "Methods")
# p_ARI
ggsave( plot = p_ARI,
        filename = "figures/Figure12.pdf",
        height = 6, width = 8)

####### Figure 13 #######

df_power %>%
  filter(n %in% c(100)) %>%
  filter(truesignal != "0") %>%
  filter(!(method_inference %in% c("tn_test", "Naif - Z.test", "Hivert - direct_test - estimated"))) %>%
  ggplot(aes(x = a, y = power, color = code)) +
  geom_point() + geom_line() + geom_abline(intercept = 0.05, slope = 0) +
  facet_wrap( vars(TV),
              labeller = label_parsed
  )+
  scale_color_manual(values = color_palette)  +
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(y = "Statistical power (through 500 experiments)", color = "Methods") -> p_power
# p_power
ggsave(plot = p_power,
       filename = "figures/Figure13.pdf",
       width = 8, height = 4)
