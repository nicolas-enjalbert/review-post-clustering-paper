pipeline_t.test <- function(X, clustering, comparison, MatMean,
                            binary_comparison,
                            schema_comparison = c("1vs1", "1vsothers")) {
  p <- dim(X)[2]
  nb_comparison <- dim(comparison)[1]

  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )

  matrix_pvalues <- matrix(NA, nrow = p, ncol = nb_comparison)
  vector_time_pvalue <- c()
  time_inference_start <- Sys.time()
  for (k in 1:nb_comparison) {
    if (all(comparison[k, ] %in% clustering)) {
      for (j in 1:p) {
        time_pvalue_start <- Sys.time()
        Xk <- X[which(clustering == comparison[k, 1]), j]
        Xothers <- X[which(clustering %in% comparison[k, -1]), j]
        if (min(length(Xk), length(Xothers)) >= 2) {
          matrix_pvalues[j, k] <- t.test(Xk, Xothers)$p.value
        }
        time_pvalue_stop <- Sys.time()
        vector_time_pvalue <- c(
          vector_time_pvalue,
          time_pvalue_stop - time_pvalue_start
        )
      }
    }
  }
  time_inference_stop <- Sys.time()

  c(
    list(
      matrix_pvalues = matrix_pvalues,
      time_mean_pvalue = mean(vector_time_pvalue),
      time_inference = time_inference_stop - time_inference_start,
      comparison = comparison,
      parameters = c(test = "t.test")
    ),
    clustering_comparison_truth
  )
}


pipeline_Z.test <- function(X, Sigma, clustering, comparison, MatMean,
                            binary_comparison,
                            schema_comparison = c("1vs1", "1vsothers")) {
  p <- dim(X)[2]
  nb_comparison <- dim(comparison)[1]

  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )

  matrix_pvalues <- matrix(NA, nrow = p, ncol = nb_comparison)
  vector_time_pvalue <- c()
  time_inference_start <- Sys.time()
  for (k in 1:nb_comparison) {
    if (all(comparison[k, ] %in% clustering)) {
      eta <- eta_comparison_mean(clustering,
        k1 = comparison[k, 1],
        k2 = comparison[k, -1]
      )
      for (j in 1:p) {
        time_pvalue_start <- Sys.time()
        Xj <- X[, j]
        sig <- Sigma[j, j]
        matrix_pvalues[j, k] <- univariate_Z_test(Xj,
          eta = eta,
          sig = sqrt(sig)
        )$pval
        time_pvalue_stop <- Sys.time()
        vector_time_pvalue <- c(
          vector_time_pvalue,
          time_pvalue_stop - time_pvalue_start
        )
      }
    }
  }
  time_inference_stop <- Sys.time()
  c(
    list(
      matrix_pvalues = matrix_pvalues,
      time_mean_pvalue = mean(vector_time_pvalue),
      time_inference = time_inference_stop - time_inference_start,
      comparison = comparison,
      parameters = c(test = "Z.test")
    ),
    clustering_comparison_truth
  )
}


pipeline_hivert <- function(X, clustering, sigma,
                            comparison, ndraws = 1000,
                            type_test = c(
                              "direct_test", "merge_test",
                              "multimod_test"
                            ),
                            clustering_method = c("kmeans", "HAC"),
                            MatMean, binary_comparison,
                            schema_comparison = c("1vs1", "1vsothers")) {
  type_test <- match.arg(type_test)

  n <- dim(X)[1]
  p <- dim(X)[2]
  K_fix <- length(unique(clustering))
  nb_comparison <- dim(comparison)[1]


  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )


  matrix_pvalues <- matrix(NA, nrow = p, ncol = nb_comparison)
  vector_time_pvalue <- c()
  time_inference_start <- Sys.time()
  for (k in 1:nb_comparison) {
    for (j in 1:p) {
      time_pvalue_start <- Sys.time()

      if (clustering_method == "kmeans") {
        cl_fun <- function(X) {
          n <- dim(X)[1]
          clustering_a <- pipeline_kmeans(X, K_fix)$clustering
          clustering <- rep(3, n)
          clustering[which(clustering_a %in% comparison[k, 1])] <- 1
          clustering[which(clustering_a %in% comparison[k, -1])] <- 2
          return(clustering)
        }
      } else if (clustering_method == "HAC") {
        cl_fun <- function(X) {
          n <- dim(X)[1]
          clustering_a <- pipeline_HAC(X, K_fix)$clustering
          clustering <- rep(3, n)
          clustering[which(clustering_a %in% comparison[k, 1])] <- 1
          clustering[which(clustering_a %in% comparison[k, -1])] <- 2
          return(clustering)
        }
      }

      clust_arrange <- rep(3, n)
      clust_arrange[which(clustering %in% comparison[k, 1])] <- 1
      clust_arrange[which(clustering %in% comparison[k, -1])] <- 2

      if (type_test == "direct_test") {
        matrix_pvalues[j, k] <- test_selective_inference(
          X = X, k1 = 1,
          k2 = 2, g = j,
          cl_fun = cl_fun,
          cl = clust_arrange,
          ndraws = ndraws,
          sig = sigma
        )$pval
      } else if (type_test == "merge_test") {
        matrix_pvalues[j, k] <- merge_selective_inference(
          X = X, k1 = 1,
          k2 = 2, g = j,
          cl_fun = cl_fun,
          cl = clust_arrange,
          ndraws = ndraws
        )$pval
      } else if (type_test == "multimod_test") {
        matrix_pvalues[j, k] <- test_multimod(
          X = X, k1 = 1,
          k2 = 2, g = j,
          cl = clust_arrange
        )$pval
      }



      time_pvalue_stop <- Sys.time()
      vector_time_pvalue <- c(
        vector_time_pvalue,
        time_pvalue_stop - time_pvalue_start
      )
    }
  }
  time_inference_stop <- Sys.time()

  c(
    list(
      matrix_pvalues = matrix_pvalues,
      time_mean_pvalue = mean(vector_time_pvalue),
      time_inference = time_inference_stop - time_inference_start,
      comparison = comparison,
      parameters = c(test = type_test)
    ),
    clustering_comparison_truth
  )
}


pipeline_data_thinning <- function(X, K, a, Sigma, epsilon = 0.5,
                                   comparison,
                                   clustering_method = c("kmeans", "HAC"),
                                   test = c("t.test", "wilcox.test"),
                                   MatMean, binary_comparison,
                                   schema_comparison = c("1vs1", "1vsothers")) {
  test <- match.arg(test)
  thinning <- datathin_gaussian_multivariate(X, epsilon, Sigma)

  Xclustering <- thinning$Xtr
  Xinference <- thinning$Xte

  if (clustering_method == "kmeans") {
    res_clustering <- pipeline_kmeans(Xclustering, K)
    clustering <- res_clustering$clustering
    time_clustering <- res_clustering$time_clustering
  } else if (clustering_method == "HAC") {
    res_clustering <- pipeline_HAC(Xclustering, K)
    clustering <- res_clustering$clustering
    time_clustering <- res_clustering$time_clustering
  }

  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )

  if (test == "t.test") {
    res_inference <- pipeline_t.test(
      X = Xinference, clustering = clustering,
      comparison = comparison,
      MatMean = MatMean,
      binary_comparison = binary_comparison,
      schema_comparison = schema_comparison
    )
  } else if (test == "Z.test") {
    res_inference <- pipeline_Z.test(
      X = Xinference, Sigma = (1 - epsilon) * Sigma,
      clustering = clustering,
      comparison = comparison,
      MatMean = MatMean,
      binary_comparison = binary_comparison,
      schema_comparison = schema_comparison
    )
  }
  ## Add here other pipeline of statistic test

  return(c(res_clustering, res_inference, clustering_comparison_truth))
}


pipeline_chen_marginal <- function(X, K, a, Sigma,
                                   comparison,
                                   clustering,
                                   clustering_method = c("kmeans", "HAC"),
                                   MatMean, binary_comparison,
                                   schema_comparison = c("1vs1", "1vsothers")) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  K_fix <- length(unique(clustering))
  nb_comparison <- dim(comparison)[1]


  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )

  if (clustering_method == "HAC") {
    hcl <- hclust(dist(X, method = "euclidean")^2, method = "ward.D")
  }

  time_inference_start <- Sys.time()
  matrix_pvalues <- matrix(NA, nrow = p, ncol = nb_comparison)
  vector_time_pvalue <- c()
  time_inference_start <- Sys.time()
  for (k in 1:nb_comparison) {
    for (j in 1:p) {
      time_pvalue_start <- Sys.time()

      if (clustering_method == "kmeans") {
        matrix_pvalues[j, k] <- kmeans_inference_1f(
          X = X, k = K_fix,
          cluster_1 = comparison[k, 1],
          cluster_2 = comparison[k, 2],
          feat = j,
          iso = FALSE,
          covMat = Sigma
        )$pval
      } else if (clustering_method == "HAC") {
        matrix_pvalues[j, k] <- test_hier_clusters_exact_1f(
          X = X,
          link = "ward.D",
          hcl = hcl,
          K = K_fix,
          k1 = comparison[k, 1],
          k2 = comparison[k, 2],
          feat = j,
          indpt = FALSE,
          covMat = Sigma
        )$pval
      }



      time_pvalue_stop <- Sys.time()
      vector_time_pvalue <- c(
        vector_time_pvalue,
        time_pvalue_stop - time_pvalue_start
      )
    }
  }
  time_inference_stop <- Sys.time()

  return(c(
    list(
      matrix_pvalues = matrix_pvalues,
      time_mean_pvalue = mean(vector_time_pvalue),
      time_inference = time_inference_stop - time_inference_start,
      comparison = comparison,
      parameters = c(test = "t.test")
    ),
    clustering_comparison_truth
  ))
}


pipeline_poclin <- function(X, K, Sigma, aggregation = c("rescale", "walesiak"),
                            comparison,
                            MatMean, binary_comparison,
                            schema_comparison = c("1vs1", "1vsothers")) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  nb_comparison <- dim(comparison)[1]

  time_clustering_start <- Sys.time()
  X_null <- SimuMultiCorr(n = n, MatMean = 0, Sigma = diag(x = 1, 1000))

  ## Choose value of lambda
  lambdas <- get_lambda_max(X_null)
  lambda <- min(lambdas) - sd(lambdas)

  trees <- convex_clustering(X)
  clust <- cutree_ordinal(trees = trees, h = lambda)

  if (aggregation == "rescale") {
    ## Aggregation on rescaled data
    clust_rescaled <- sweep(clust - 1, 2, matrixStats::colMaxs(clust - 1), "/")
    tree <- hclust(dist(clust_rescaled, "euclidean"), "ward.D2")
  } else if (aggregation == "walesiak") {
    ## Aggregation with Walesiak metric
    dist <- GDM2(clust)
    tree <- hclust(dist, "ward.D2")
  }
  clustering <- cutree(tree, k = K)
  time_clustering_stop <- Sys.time()


  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )

  Gamma <- Sigma %x% Diagonal(n, rep(1, n))
  matrix_pvalues <- matrix(NA, nrow = p, ncol = nb_comparison)
  time_inference_start <- Sys.time()
  for (k in 1:nb_comparison) {
    eta <- rep(0, n)
    eta[clustering == comparison[k, 1]] <- 1 / sum(clustering == comparison[k, 1])
    eta[clustering %in% comparison[k, -1]] <- -1 / sum(clustering %in%
      comparison[k, -1])
    for (j in 1:p) {
      ind <- 1:p == j
      kappa <- kronecker(eta, t(ind))
      matrix_pvalues[j, k] <- get_cond_p_value(
        X = X, kappa = kappa,
        clust = clust,
        lambda = lambda,
        Gamma = as.matrix(Gamma)
      )
    }
  }
  time_inference_stop <- Sys.time()

  c(
    list(
      matrix_pvalues = matrix_pvalues,
      time_mean_pvalue = (time_inference_stop - time_inference_start) / (K * p),
      time_inference = time_inference_stop - time_inference_start,
      comparison = comparison,
      parameters = c(
        test = "Poclin", method = "poclin",
        clustering_method = aggregation
      ),
      clustering = clustering,
      time_clustering = time_clustering_stop - time_clustering_start
    ),
    clustering_comparison_truth
  )
}

pipeline_tn_test <- function(X, K, comparison,
                             clustering_method = c("kmeans", "HAC", "GMM"),
                             MatMean, binary_comparison,
                             schema_comparison = c("1vs1", "1vsothers")) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  nb_comparison <- dim(comparison)[1]

  ## Splitting
  inds1 <- sort(sample(1:n, size = n %/% 2, replace = FALSE))
  Xclustering <- X[inds1, ]
  Xinference <- X[-inds1, ]

  samp_labels <- rep(0, n)
  samp_labels[inds1] <- 1
  samp_labels <- as.factor(samp_labels)

  ## Clustering on Xclustering
  if (clustering_method == "kmeans") {
    res_clustering <- pipeline_kmeans(Xclustering, K)
    clustering <- res_clustering$clustering
    time_clustering <- res_clustering$time_clustering
  } else if (clustering_method == "HAC") {
    res_clustering <- pipeline_HAC(Xclustering, K)
    clustering <- res_clustering$clustering
    time_clustering <- res_clustering$time_clustering
  }

  ## Inference on
  truncated_normal <- import("truncated_normal")
  tn_test <- truncated_normal$truncated_normal
  matrix_pvalues <- matrix(NA, nrow = p, ncol = nb_comparison)
  time_inference_start <- Sys.time()
  if (dim(comparison)[2] > 2) {
    k <- 1
    for (k in comparison[, 1]) {
      contrast <- rep(0, length(clustering))
      contrast[which(clustering == k)] <- 1
      SVC_k <- tn_test$SVC(kernel = "linear", C = 100)
      SVC_k$fit(Xclustering, contrast)
      labelsInference_k <- SVC_k$predict(Xinference)

      y <- matrix(Xinference[which(labelsInference_k == 1), ], ncol = p)
      z <- matrix(Xinference[which(labelsInference_k == 0), ], ncol = p)
      a <- array_reshape(SVC_k$coef_, dim = -1)
      b <- SVC_k$intercept_
      pvals <- tn_test$tn_test(
        y = y, z = z, a = a, b = b,
        learning_rate = 1, eps = 1e-2,
        verbose = FALSE, return_likelihood = FALSE,
        num_iters = as.integer(10000), num_cores = NULL
      )
      matrix_pvalues[, k] <- pvals
    }
  } else if (dim(comparison)[2] == 2) {
    #### Fit hyperplanes using X1
    SVC <- tn_test$SVC(kernel = "linear", C = 100)
    SVC$fit(Xclustering, clustering)

    ### General labels using hyperplanes on Xinference
    labelsInference <- SVC$predict(Xinference)
    for (k in 1:nb_comparison) {
      y <- matrix(Xinference[which(labelsInference == comparison[k, 1]), ],
        ncol = p
      )
      z <- matrix(Xinference[which(labelsInference == comparison[k, 2]), ],
        ncol = p
      )
      a <- array_reshape(SVC$coef_[k, ], dim = -1)
      b <- SVC$intercept_[k]
      pvals <- tn_test$tn_test(
        y = y, z = z, a = a, b = b,
        learning_rate = 1, eps = 1e-2,
        verbose = FALSE, return_likelihood = FALSE,
        num_iters = as.integer(10000), num_cores = NULL
      )
      matrix_pvalues[, k] <- pvals
    }
  }
  time_inference_stop <- Sys.time()

  clustering_final <- rep(NA, n)
  clustering_final[inds1] <- clustering
  clustering_final[-inds1] <- labelsInference

  clustering_comparison_truth <- truth_comparison_clusters(
    clustering = clustering_final,
    MatMean = MatMean,
    binary_comparison = binary_comparison,
    schema_comparison = schema_comparison
  )

  c(
    list(
      matrix_pvalues = matrix_pvalues,
      time_mean_pvalue = (time_inference_stop - time_inference_start) /
        (nb_comparison * p),
      time_inference = time_inference_stop - time_inference_start,
      comparison = comparison,
      parameters = c(
        test = "Tn_test", method = "Tn_test",
        clustering_method = clustering_method
      ),
      clustering = clustering_final,
      time_clustering = time_clustering
    ),
    clustering_comparison_truth
  )
}
