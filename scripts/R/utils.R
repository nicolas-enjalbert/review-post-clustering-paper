library(tidyverse)
## Simulation functions ##

signal_3clusters_eq <- function(n, a = 10, p = 3,
                                prob_clusters = c(1 / 3, 1 / 3, 1 / 3),
                                verbose = TRUE) {
  K <- 3
  truthF <- rep(0, p)
  truthS <- sort(c(1:K, sample(1:K,
                  size = n - K, replace = TRUE,
                  prob = prob_clusters
  )))

  X <- matrix(0, nrow = n, ncol = p)
  signal <- matrix(c(
    -a / 2, rep(0, p - 1),
    a / 2, rep(0, p - 1),
    0, sqrt(3) * a / 2, rep(0, p - 2)
  ), nrow = 3, byrow = TRUE)
  for (k in 1:K) {
    indS <- which(truthS == k)
    X[indS, ] <- rep(signal[k, ], each = length(indS))
  }
  if (verbose) plot(X)

  return(list(X = X, truthF = truthF, truthS = truthS))
}

signal_2clusters <- function(n, p, a = 10, prob_clusters = c(0.5, 0.5),
                             nb_variable_signal = p, verbose = TRUE) {
  K <- 2
  truthF <- c(rep(1, nb_variable_signal), rep(0, p - nb_variable_signal))
  truthS <- sort(c(1:K, sample(1:K,
                  size = n - 2, replace = TRUE,
                  prob = prob_clusters
  )))

  X <- matrix(0, nrow = n, ncol = p)
  signal <- matrix(c(
    rep(0, p),
    rep(a, nb_variable_signal),
    rep(0, p - nb_variable_signal)
  ), nrow = 2, byrow = TRUE)
  for (k in 1:K) {
    indS <- which(truthS == k)
    X[indS, ] <- rep(signal[k, ], each = length(indS))
  }
  if (verbose) plot(X[, 5:6], col = truthS)
  return(list(X = X, truthF = truthF, truthS = truthS))
}

signal_2clusters_indep_p <- function(n, p, a = 10, prob_clusters = c(0.5, 0.5),
                                     nb_variable_signal = p, verbose = TRUE) {
  K <- 2
  truthF <- c(rep(1, nb_variable_signal), rep(0, p - nb_variable_signal))
  truthS <- sort(c(1:K, sample(1:K,
                               size = n - 2, replace = TRUE,
                               prob = prob_clusters
  )))

  X <- matrix(0, nrow = n, ncol = p)
  signal <- matrix(c(
    rep(0, p),
    rep(a / sqrt(nb_variable_signal), nb_variable_signal),
    rep(0, p - nb_variable_signal)
  ), nrow = 2, byrow = TRUE)
  for (k in 1:K) {
    indS <- which(truthS == k)
    X[indS, ] <- rep(signal[k, ], each = length(indS))
  }
  if (verbose) plot(X[, 1:2], col = truthS)
  return(list(X = X, truthF = truthF, truthS = truthS))
}


SimuMultiCorr <- function(n, MatMean, Sigma) {
  library(mvtnorm)
  p <- nrow(Sigma)
  X <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
  X <- X + MatMean
  return(X)
}


#### Data thinning function ####
datathin_gaussian_multivariate <- function(X, epsilon, Sigma) {
  i <- 1
  Xtr <- matrix(ncol = dim(X)[2], nrow = dim(X)[1])
  for (i in 1:dim(X)[1]) {
    Xtr[i, ] <- MASS::mvrnorm(
      n = 1, mu = epsilon * X[i, ],
      Sigma = epsilon * (1 - epsilon) * Sigma
    )
  }
  Xte <- X - Xtr
  Xdt <- list(Xtr = Xtr, Xte = Xte)
  return(Xdt)
}



eta_comparison_mean <- function(clustering, k1, k2) {
  eta <- (clustering == k1) * 1 / length(which(clustering == k1)) -
    (clustering == k2) * 1 / length(which(clustering == k2))
  return(eta)
}

norm_vec <- function(x) sqrt(sum(x^2))

multivariate_Z_test <- function(X, eta, sig) {
  p <- ncol(X)

  diff_means <- eta %*% X
  stat <- norm_vec(diff_means)
  squared_norm_eta <- norm_vec(eta)^2
  scale_factor <- squared_norm_eta * sig^2

  pval <- 1 - pchisq((stat^2) / scale_factor, df = p)
  return(list(stat = stat, pval = pval))
}

univariate_Z_test <- function(Xj, eta, sig) {
  stat <- eta %*% Xj
  sd <- sig * norm_vec(eta)
  quant <- pnorm(stat, mean = 0, sd = sd)
  pval <- 2 * min(quant, 1 - quant)
  return(list(stat = stat, pval = pval))
}



truth_comparison_clusters <- function(clustering, MatMean, binary_comparison,
                                      schema_comparison = c("1vs1", "1vsothers")
                                      ) {
  n <- length(clustering)
  if (schema_comparison == "1vs1") {
    data.frame(clustering = clustering, n = 1, indices = 1:n) %>%
      pivot_wider(
        names_from = clustering, values_from = c(n),
        values_fill = list(n = 0)
      ) %>%
      dplyr::select(-c(indices)) %>%
      as.matrix() -> binary_n_K
    binary_p_n <- t((MatMean != 0) * 1)
    mean_p_K <- t(t(binary_n_K) %*% MatMean / colSums(binary_n_K))

    truthF_comparison <- (abs(
      mean_p_K %*% binary_comparison[colnames(mean_p_K), ]))
  } else if (schema_comparison == "1vsothers") {
    truthF_comparison <- truthF
  }
  return(list(mean_p_K = mean_p_K, binary_n_K = binary_n_K,
              truthF_comparison = truthF_comparison))
}


clean_data <- function(df) {
  df_clean <- df %>%
    mutate(method_clustering = replace(method_clustering,
                                       method_clustering == "ward", "HAC")) %>%
    mutate(method_inference_abre = recode(method_inference,
                                          "naif" = "Naive",
                                          "datathinning" = "DT",
                                          "Gao" = "Cond",
                                          "unknown_variance" = "YFB",
                                          "general_dependance" = "Cond-GD",
                                          .default = levels(method_inference)
    )) %>%
    mutate(method_clustering = recode(method_clustering,
                                      "HAC_ward.D" = "HAC"
    )) %>%
    mutate(method_clustering_abre = recode(method_clustering,
                                           "HAC" = "HAC",
                                           "HAC_ward.D" = "HAC",
                                           "kmeans" = "KM",
                                           "GMM" = "GMM",
    )) %>%
    mutate(method_clustering = recode(method_clustering,
                                      kmeans = "K-means",
                                      .default = method_clustering
    )) %>%
    mutate(type_estimation_abre = recode(type_estimation,
                                         "all" = "all",
                                         "estimation_javi" = "all",
                                         "known" = "oracle",
                                         "intra" = "intra",
                                         "unknown" = ""
    )) %>%
    # mutate(method_clustering_abre = as.factor(method_clustering_abre)) %>%
    mutate(code = paste(method_inference_abre, method_clustering_abre,
                        sep = "-")) %>%
    mutate(code = ifelse(ndraws != 0, paste(code, "IS", sep = "-"), code)) %>%
    mutate(code_inf_IS = ifelse(ndraws != 0,
                                paste(method_inference_abre, "IS", sep = "-"),
                                method_inference_abre
    )) %>%
    mutate(method_clustering = replace(method_clustering,
                                       method_clustering == "HAC",
                                       "HAC (Ward's linkage)")) %>%
    mutate(method_clustering = factor(method_clustering,
                                      levels = c("HAC (Ward's linkage)",
                                                 "K-means", "GMM")
    ))

  df_clean %>%
    dplyr::select(method_inference_abre, method_clustering_abre, ndraws,
                  code) %>%
    mutate(IS = ndraws != 0) %>%
    distinct() %>%
    mutate(color = recode(method_inference_abre,
                          "Cond" = "#7CAE00",
                          "DT" = "#F8766D",
                          "Naive" = "#C77CFF",
                          "YFB" = "#00BFC4",
                          .default = "yellow"
    )) %>%
    mutate(linetype = ifelse(IS, "dotted", "solid")) -> df_color_linetype

  color_palette <- setNames(df_color_linetype$color, df_color_linetype$code)
  linetype_palette <- setNames(df_color_linetype$linetype,
                               df_color_linetype$code)


  df_color_linetype_inf_IS <- df_clean %>%
    dplyr::select(method_inference_abre, ndraws, code_inf_IS) %>%
    mutate(IS = ndraws != 0) %>%
    distinct() %>%
    mutate(color = recode(method_inference_abre,
                          "Cond" = "#7CAE00",
                          "DT" = "#F8766D",
                          "Naive" = "#C77CFF",
                          "YFB" = "#00BFC4",
                          .default = "yellow"
    )) %>%
    mutate(linetype = ifelse(IS, "dotted", "solid"))
  color_palette_inf_IS <- setNames(
    df_color_linetype_inf_IS$color,
    df_color_linetype_inf_IS$code_inf_IS
  )
  linetype_palette_inf_IS <- setNames(
    df_color_linetype_inf_IS$linetype,
    df_color_linetype_inf_IS$code_inf_IS
  )

  return(list(
    df_clean = df_clean,
    color_palette = color_palette,
    linetype_palette = linetype_palette,
    color_palette_inf_IS = color_palette_inf_IS,
    linetype_palette_inf_IS = linetype_palette_inf_IS
  ))
}

clean_data_estimation <- function(df) {
  df %>%
    filter(ndraws == 0) %>%
    mutate(method_inference_abre = recode(method_inference,
                                          "naif" = "Naive",
                                          "datathinning" = "DT",
                                          "Gao" = "Cond",
                                          "unknown_variance" = "YFB",
                                          "general_dependance" = "Cond-GD",
                                          .default = levels(method_inference)
    )) %>%
    # mutate(method_inference_abre = as.factor(method_inference_abre)) %>%
    mutate(method_clustering_abre = recode(method_clustering,
                                           "HAC" = "HAC",
                                           "kmeans" = "KM",
                                           "GMM" = "GMM"
    )) %>%
    mutate(type_estimation_abre = recode(type_estimation,
                                         "all" = "all",
                                         "estimation_javi" = "all",
                                         "exact" = "oracle",
                                         "known" = "oracle",
                                         "intra" = "intra",
                                         "unknown" = ""
    )) %>%
    # mutate(method_clustering_abre = as.factor(method_clustering_abre)) %>%
    mutate(code = ifelse(method_inference_abre == "YFB", method_inference_abre,
                         paste(method_inference_abre, type_estimation_abre,
                               sep = "-"))) %>%
    mutate(code = ifelse(ndraws != 0, paste(code, "IS", sep = "-"), code)) %>%
    mutate(code_inf_IS = ifelse(ndraws != 0,
                                paste(method_inference_abre, "IS", sep = "-"),
                                method_inference_abre
    )) %>%
    mutate(type_estimation = factor(type_estimation_abre,
                                    levels = c("oracle", "all", "intra", "")
    )) %>%
    filter(code != "Cond-GD-oracle") -> df_clean

  df_clean %>%
    dplyr::select(
      method_inference_abre, method_clustering_abre,
      type_estimation_abre, ndraws, code
    ) %>%
    mutate(IS = ndraws != 0) %>%
    distinct() %>%
    mutate(color = recode(method_inference_abre,
                          "Cond" = "#7CAE00",
                          "DT" = "#F8766D",
                          "Naive" = "#C77CFF",
                          "YFB" = "#00BFC4",
                          "Cond-GD" = "orange",
                          .default = "yellow"
    )) %>%
    mutate(linetype = recode(type_estimation_abre,
                             "oracle" = "solid",
                             "intra" = "twodash",
                             "all" = "dotted",
                             "IS" = "dotted",
                             .default = "solid"
    )) -> df_color_linetype

  color_palette <- setNames(df_color_linetype$color, df_color_linetype$code)
  linetype_palette <- setNames(df_color_linetype$linetype,
                               df_color_linetype$code)

  return(list(
    df_clean = df_clean,
    color_palette = color_palette,
    linetype_palette = linetype_palette
  ))
}


clean_data_marginal <- function(df) {
  df %>%
    mutate(method_clustering = recode(method_clustering,
                                      "HAC_ward.D" = "HAC",
                                      .default = method_clustering
    )) %>%
    mutate(code = recode(method_inference,
                         "Naif - Z.test" = "Naive",
                         "tn_test" = "TN",
                         "Data thinning - Z.test" = "DT",
                         "Chen-marginal" = "CG",
                         "Hivert - direct_test" = "H",
                         "Hivert - merge_test" = "H-merge",
                         "poclin" = "poclin",
                         "Hivert - multimod_test" = "Dip",
                         .default = method_inference
    )) %>%
    group_by(
      nexp, K, n, p, ndraws, a, epsilon, method_inference, method_clustering,
      exact, sigma, rho, variables
    ) %>%
    mutate(true_C1_vs_C2 = case_when(
      variables == 1 ~ max(truehypothesis, na.rm = TRUE),
      variables == 2 ~ min(truehypothesis, na.rm = TRUE),
      .default = min(truehypothesis, na.rm = TRUE)
    )) %>%
    ungroup() %>%
    mutate(H0 = case_when(
      variables == 3 ~ "H0 j = 3",
      variables == 1 ~ "H1 j = 1",
      variables == 2 & truehypothesis == true_C1_vs_C2 ~ "H0 j = 2",
      variables == 2 & truehypothesis != true_C1_vs_C2 ~ "H1 j = 2"
    )) %>%
    mutate(truesignal = case_when(
      variables == 3 ~ "0",
      variables == 1 & truehypothesis == true_C1_vs_C2 ~ "a",
      variables == 1 & truehypothesis != true_C1_vs_C2 ~ "a/2",
      variables == 2 & truehypothesis == true_C1_vs_C2 ~ "0",
      variables == 2 & truehypothesis != true_C1_vs_C2 ~ "$\\sqrt{3}a/2$"
    )) %>%
    mutate(truesignal_numeric = case_when(
      variables == 3 ~ 0,
      variables == 1 & truehypothesis == true_C1_vs_C2 ~ a,
      variables == 1 & truehypothesis != true_C1_vs_C2 ~ a / 2,
      variables == 2 & truehypothesis == true_C1_vs_C2 ~ 0,
      variables == 2 & truehypothesis != true_C1_vs_C2 ~ sqrt(3) * a / 2
    )) %>%
    mutate(TV = case_when(
      truesignal == "a" & variables == 1 ~ 1,
      truesignal == "a/2" & variables == 1 ~ 2,
      truesignal == "$\\sqrt{3}a/2$" & variables == 2 ~ 3
    )) %>%
    mutate(TV = factor(TV,
                       labels = c(
                         "1" = expression("signal: " * a * ",    " * j == 1),
                         "2" = expression(paste("signal: ", frac(a, 2), ",    ",
                                                j == 1)),
                         "3" = expression(paste("signal: ",
                                                frac(sqrt(3) * a, 2), ",    ",
                                                j == 2))
                       )
    )) -> df_clean

  df_clean %>%
    dplyr::select(code) %>%
    distinct() %>%
    arrange(code) -> df_color_linetype
  color <- c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3",
             "#619CFF", "#DB72FB", "#FF61C3")
  color_palette <- setNames(color, df_color_linetype$code)

  return(list(
    df_clean = df_clean,
    color_palette = color_palette
  ))
}

estimation_cov <- function(X, cl, rho) {
  K <- length(unique(cl))
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (rho == 0) {
    Sigma_all <- sum(rowMeans((X - colMeans(X))^2)) / (n - 1) * diag(p)

    vec_sigma_k <- c()
    for (k in 1:K) {
      Xk <- X[cl == k, ]
      nk <- length(which(cl == k))
      sigma_k <- sum(rowMeans(
        (scale(Xk, center = TRUE, scale = FALSE))^2
      )) / (nk - 1)
      vec_sigma_k[k] <- sigma_k * nk
    }
    sigma_intra <- sum(vec_sigma_k) / n
    Sigma_intra <- sigma_intra * diag(p)
  } else {
    Sigma_all <- cov(X)
    covmatrix <- matrix(0, ncol = p, nrow = p)
    for (k in 1:K) {
      nk <- length(which(cl == k))
      covmatrix <- covmatrix + cov(X[cl == k, ]) * nk
    }
    Sigma_intra <- covmatrix / n
  }
  return(list(Sigma_all = Sigma_all, Sigma_intra = Sigma_intra))
}
