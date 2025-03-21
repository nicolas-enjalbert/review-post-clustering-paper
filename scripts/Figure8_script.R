library(tidyverse)
source("scripts/R/script_simulation_multivariate-methods.R")

setting <- 1
path <- sprintf("results/simulation_multivariate_setting%s/Figure8/", setting)

####### Compute experiments #######
comparison_multivariate_inference_methods(
  methods = c("Cond", "YFB"),
  clusterings = c("HAC"),
  nb_experiments = 500,
  setting = setting,
  ns = c(500),
  ps = c(2),
  as = c(0:10),
  Ndraws = c(0, 1000),
  epsilon_dts = c(0.7),
  sigmas = c(1),
  rhos = c(0.3, 0.5),
  type_estimation = c("known", "all", "intra"),
  path = path,
  nb_workers = min(40, availableCores() / 4)
)

####### Figures #######

pattern <- "sim_multi_set1_K=2_a=(.*)_n=(.*)_p=(.*)_nb_experiment=(.*)_sigma=(.*)_rho=(.*).RDS"
filenames <- list.files(path, pattern = pattern)

df_end <- data.frame()
for (filename in filenames) {
  df <- readRDS(paste(path, filename, sep = "/"))
  df_end <- rbind(df_end, df)
}

data_cleaned <- clean_data_estimation(df_end)
df_clean <- data_cleaned$df_clean
color_palette <- data_cleaned$color_palette
linetype_palette <- data_cleaned$linetype_palette

####### Figure 8 #######


p_ecdf_dep <- df_clean %>%
  filter(a == 0) %>%
  filter(code != c("Cond-all")) %>%
  ggplot(aes(x = pvaleur, color = code, linetype = code)) +
  stat_ecdf(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(rho ~ p,
    labeller = label_bquote(
      rows = rho ~ "=" ~ .(rho),
      cols = m ~ "=" ~ .(p)
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  labs(x = "p-values", color = "Methods", linetype = "Methods")

p_power_dep <- df_clean %>%
  filter(p %in% c(2)) %>%
  filter(n == 500) %>%
  filter(code %in% c("Cond-GD-all", "Cond-oracle")) %>%
  mutate(reject = pvaleur <= 0.05) %>%
  group_by(
    n, p, a, ndraws, method_inference, method_clustering, exact, sigma,
    rho, type_estimation, code
  ) %>%
  summarise(
    nb_exp = n(),
    power = mean(reject, na.rm = TRUE),
    mean_time = mean(time, na.rm = TRUE),
    mean_ARI = mean(ARI, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = a, y = power, color = code, linetype = code)) +
  geom_point() +
  geom_line(linewidth = 1.1) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text.x = element_blank()
  ) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  labs(color = "Methods", linetype = "Methods", y = "Statistical power") +
  facet_grid(rho ~ p,
    labeller = label_bquote(
      rows = rho ~ "=" ~ .(rho),
      cols = m ~ "=" ~ .(p)
    )
  )

figure8 <- ggpubr::ggarrange(p_ecdf_dep, p_power_dep,
  labels = c("A", "B"),
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom"
)

ggsave(
  plot = figure8,
  "figures/Figure8.pdf",
  width = 8, height = 4
)
