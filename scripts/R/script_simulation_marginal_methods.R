library(future)
library(future.apply)
library(tidyverse)
library(Rmixmod)
library(VALIDICLUST)
library(CADET)
library(poclin)
library(Matrix)
library(clusterSim)
library(reticulate)
library(future.callr)
source("scripts/R/utils.R")
source("scripts/R/pipeline_clustering.R")
source("scripts/R/pipeline_inference.R")

#' Compute marginal post-clustering inference methods
#'
#' @param methods vector of strings, the evaluated marginal post-clustering inference methods
#' @param clusterings vector of strings, the used clustering methods
#' @param nb_experiments int, number of experiments
#' @param setting int, setting number used
#' @param ns vector of integers, numbers of samples used in simulations
#' @param ps vector of integers, numbers of variables used in simulations
#' @param as vector of integers, signals used in simulations
#' @param Ndraws vector of integers, numbers of Monte Carlo Importance Sampling draws used in Hivert's methods
#' @param epsilon_dts vector of float, values of epsilon between 0 and 1 used in Data thinning methods.
#' @param sigmas vector of float, values of the variance
#' @param rhos vector of float, values of the covariance
#' @param path string, path to save the results
#' @param nb_workers number of workers for the future package to parallelize the computation
#'
#' @return
comparison_univariate_inference_methods <- function(
    methods = c("Z.test", "data.thinning", "Hivert", "Chen", "tn.test", "poclin"),
    clusterings = c("HAC", "kmeans", "GMM"),
    nb_experiments = 500, setting = c(1,2),
    ns = 100,
    ps = 3,
    as = c(0,10),
    Ndraws = c(1000),
    epsilon_dts = c(0.1, 0.3, 0.5, 0.7, 0.9),
    sigmas = c(1),
    rhos = c(0),
    path =sprintf("results/simulation_marginal_setting%s", setting),
    nb_workers = min(40, availableCores()/4)
){


  plan(multisession, workers = nb_workers)
  set.seed(130124)


  dir.create(path, showWarnings = FALSE, recursive = TRUE)


  configs <- expand.grid(ns = ns,
                         as = as,
                         ps = ps,
                         sigmas = sigmas,
                         rhos = rhos)
  seq_configs <- 1:nrow(configs)





  schema_comparison = "1vs1"
  if (setting == 1){
    K = 2
  } else if (setting == 2){
    K = 3
  }



  cc = 1
  for (cc in seq_configs) {
    config <- configs[cc,]
    a <- config[["as"]]
    n <- config[["ns"]]
    p <- config[["ps"]]
    sigma <- config[["sigmas"]]
    rho <- config[["rhos"]]

    print(sprintf("a = %s ; n = %s ; p = %s ; sigma = %s ; rho  = %s",
                  a, n, p, sigma, rho))
    {comparison <- t(combn(x = 1:K, 2))

      comparison %>% data.frame() %>%
        mutate(n = 1) %>%
        rownames_to_column(var = "indices") %>%
        pivot_wider(names_from = c(X1), values_from = c(n),
                    values_fill = list(n=0)) %>%
        mutate(n = -1) %>%
        pivot_wider(names_from = c(X2), values_from = c(n),
                    values_fill = list(n=0), names_repair = "minimal") -> df1

      binary_comparison <- rowsum(t(as.data.frame(df1[,-1])),
                                  group = colnames(as.data.frame(df1[,-1])),
                                  na.rm = T)
    }

    nexp <- 1
    res <- future.apply::future_lapply(1:nb_experiments, future.seed = TRUE, FUN = function(nexp){
      # res <- lapply(1:nb_experiments, FUN = function(nexp){
      print(nexp)


      {

        if (setting == 1){
          signal <- signal_2clusters(n = n, p = p, a = a, prob_clusters = c(0.5, 0.5),
                                     nb_variable_signal = p/2, verbose = FALSE)
          K = 2
          Sigma <- toeplitz(rho^(0:(p-1))) * sigma
          SigmaInv <- solve(Sigma)
        } else if (setting == 2){
          signal <- signal_3clusters_eq(n, a = a, p = p,
                                        prob_clusters = c(1/3,1/3,1/3),
                                        verbose = FALSE)
          K = 3
          sigmas <- rep(sigma,p)
          Sigma <- diag(sigmas)
          Sigma[1,3] <- Sigma[3,1] <- rho
        }

        MatMean <- signal$X
        truthS <- signal$truthS
        truthF <- signal$truthF

      }

      X <- SimuMultiCorr(n = n, MatMean = MatMean, Sigma = Sigma)

      clustering_kmeans <- pipeline_kmeans(X, K)
      clustering_HAC <- pipeline_HAC(X, K)
      clustering_GMM <- pipeline_GMM(X, K)

      df <- data.frame()


      ############# NAIVE METHOD #################





      # Z test + HAC
      if("Z.test" %in% methods){
        for(cl_meth in clusterings){
          if(cl_meth == "kmeans"){
            clustering <- clustering_kmeans$clustering
          } else if (cl_meth == "HAC"){
            clustering <- clustering_HAC$clustering
          } else if (cl_meth == "GMM"){
            clustering <- clustering_GMM$clustering
          }
        Z_test_HAC <- pipeline_Z.test(X = X, Sigma = Sigma,
                                      clustering = clustering,
                                      comparison = comparison,
                                      MatMean = MatMean,
                                      binary_comparison= binary_comparison,
                                      schema_comparison = schema_comparison)
        ARI <- mclust::adjustedRandIndex(truthS, clustering)

        results <- Z_test_HAC
        matrix_pvalues <- results$matrix_pvalues %>% data.frame()
        colnames(matrix_pvalues)  <- paste("X", 1:ncol(matrix_pvalues), sep = "")
        truthF_comparison <- results$truthF_comparison %>% data.frame()
        colnames(truthF_comparison) <- paste("X", 1:ncol(truthF_comparison), sep = "")
        comparison <- results$comparison
        time_ = results$time_mean_pvalue

        gene_truehypo <- ((truthF_comparison)) %>%
          data.frame() %>%
          rowid_to_column(var = "variables") %>%
          pivot_longer(cols = c(-variables), names_to = "cluster", values_to = "truehypothesis") %>%
          left_join(
            comparison %>%
              data.frame() %>%
              rowid_to_column(var = "cluster") %>%
              mutate(cluster = paste("X", cluster, sep = "")) %>%
              mutate(comparison = paste(X1, "vs", X2)) %>%
              dplyr::select(-c(X1, X2)),
            by = c("cluster")
          ) %>% mutate(variables = as.character(variables))

        pvalues_ext <- matrix_pvalues %>%
          data.frame() %>%
          rownames_to_column(var = "variables") %>%
          pivot_longer(cols = c(-variables),
                       names_to = "cluster",
                       values_to = "pvalues") |>
          left_join(gene_truehypo ,by = c("variables", "cluster"))


        df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                   ndraws = 0, a = a, epsilon = 0,
                                   method_inference = "Naif - Z.test",
                                   method_clustering = cl_meth,
                                   exact = TRUE,
                                   sigma = sigma,  rho = rho,
                                   pvalues_ext,
                                   ARI = ARI, time = as.numeric(time_ )))
        rm(Z_test_HAC)
        rm(gene_truehypo)
        rm(pvalues_ext)
        }
      }


      ############## DATA THNNING ################
      TESTS_prime <- c("t.test")

      if("data.thinning" %in% methods){
        for(cl_meth in clusterings){
          for (test in TESTS_prime){
            for(epsilon in epsilon_dts){
              data_thinning <- pipeline_data_thinning(X = X, K = K, a = a,
                                                      Sigma = Sigma, epsilon = epsilon,
                                                      comparison = comparison,
                                                      clustering_method = cl_meth,
                                                      test = test, MatMean = MatMean,
                                                      binary_comparison= binary_comparison,
                                                      schema_comparison = schema_comparison
              )


              ARI <- mclust::adjustedRandIndex(truthS, data_thinning$clustering)

              results <- data_thinning
              matrix_pvalues <- results$matrix_pvalues %>% data.frame()
              colnames(matrix_pvalues)  <- paste("X", 1:ncol(matrix_pvalues), sep = "")
              truthF_comparison <- results$truthF_comparison %>% data.frame()
              colnames(truthF_comparison) <- paste("X", 1:ncol(truthF_comparison), sep = "")
              comparison <- results$comparison
              time_ = results$time_mean_pvalue

              gene_truehypo <- ((truthF_comparison)) %>%
                data.frame() %>%
                rowid_to_column(var = "variables") %>%
                pivot_longer(cols = c(-variables), names_to = "cluster", values_to = "truehypothesis") %>%
                left_join(
                  comparison %>%
                    data.frame() %>%
                    rowid_to_column(var = "cluster") %>%
                    mutate(cluster = paste("X", cluster, sep = "")) %>%
                    mutate(comparison = paste(X1, "vs", X2)) %>%
                    dplyr::select(-c(X1, X2)),
                  by = c("cluster")
                ) %>% mutate(variables = as.character(variables))

              pvalues_ext <- matrix_pvalues %>%
                data.frame() %>%
                rownames_to_column(var = "variables") %>%
                pivot_longer(cols = c(-variables),
                             names_to = "cluster",
                             values_to = "pvalues") |>
                left_join(gene_truehypo ,by = c("variables", "cluster"))

              df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                         ndraws = 0, a = a, epsilon = epsilon,
                                         method_inference = paste("Data thinning", test, sep = " - "),
                                         method_clustering = cl_meth,
                                         exact = TRUE,
                                         sigma = sigma,  rho = rho,
                                         pvalues_ext,
                                         ARI = ARI, time = as.numeric(time_ )))

              rm(data_thinning)
              rm(gene_truehypo)
              rm(pvalues_ext)
            }
          }
        }
      }

      ############## HIVERT ######################
      if("Hivert" %in% methods){
        if(K == 2){
          TESTS_prime <- c( "direct_test", "multimod_test")
        } else if (K == 3) {
          TESTS_prime <- c( "direct_test", "merge_test", "multimod_test")
        }

        for(cl_meth in clusterings){
          if(cl_meth == "kmeans"){
            clustering <- clustering_kmeans
          } else if (cl_meth == "HAC"){
            clustering <- clustering_HAC
          } else if (cl_meth == "GMM"){
            clustering <- clustering_GMM
          }
          for (test in TESTS_prime){
            if(test == "multimod_test"){
              NDRAWS <- 0
            } else {
              NDRAWS <- Ndraws[Ndraws!=0]
            }

            for(ndraws in NDRAWS){
              hivert <- pipeline_hivert(X = X,
                                        clustering = clustering$clustering,
                                        sigma = sigma, ndraws = ndraws,
                                        type_test = test,
                                        comparison = comparison,
                                        clustering_method = cl_meth,
                                        MatMean = MatMean,
                                        binary_comparison= binary_comparison,
                                        schema_comparison = schema_comparison)
              ARI <- mclust::adjustedRandIndex(truthS, clustering$clustering)

              results <- hivert
              matrix_pvalues <- results$matrix_pvalues %>% data.frame()
              colnames(matrix_pvalues)  <- paste("X", 1:ncol(matrix_pvalues), sep = "")
              truthF_comparison <- results$truthF_comparison %>% data.frame()
              colnames(truthF_comparison) <- paste("X", 1:ncol(truthF_comparison), sep = "")
              comparison <- results$comparison
              time_ = results$time_mean_pvalue

              gene_truehypo <- ((truthF_comparison)) %>%
                data.frame() %>%
                rowid_to_column(var = "variables") %>%
                pivot_longer(cols = c(-variables), names_to = "cluster", values_to = "truehypothesis") %>%
                left_join(
                  comparison %>%
                    data.frame() %>%
                    rowid_to_column(var = "cluster") %>%
                    mutate(cluster = paste("X", cluster, sep = "")) %>%
                    mutate(comparison = paste(X1, "vs", X2)) %>%
                    dplyr::select(-c(X1, X2)),
                  by = c("cluster")
                ) %>% mutate(variables = as.character(variables))

              pvalues_ext <- matrix_pvalues %>%
                data.frame() %>%
                rownames_to_column(var = "variables") %>%
                pivot_longer(cols = c(-variables),
                             names_to = "cluster",
                             values_to = "pvalues") |>
                left_join(gene_truehypo ,by = c("variables", "cluster"))

              df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                         ndraws = ndraws, a = a, epsilon = 0,
                                         method_inference = paste("Hivert", test, sep = " - "),
                                         method_clustering = cl_meth,
                                         exact = TRUE,
                                         sigma = sigma,  rho = rho,
                                         pvalues_ext,
                                         ARI = ARI, time = as.numeric(time_ )))
              rm(hivert)
              rm(gene_truehypo)
              rm(pvalues_ext)
            }


          }
        }
      }

      ################ Chen Gao marginal test ###############

      if("Chen" %in% methods){
        for(cl_meth in intersect(clustering, c("HAC", "kmeans"))){
          {

            if(cl_meth == "HAC"){
              clustering <- clustering_HAC
            } else if (cl_meth == "kmeans"){
              clustering <- clustering_kmeans
            }

            chen_marginal <- pipeline_chen_marginal(X = X, Sigma = Sigma,
                                                    clustering = clustering$clustering,
                                                    comparison = comparison,
                                                    clustering_method = cl_meth,
                                                    MatMean = MatMean,
                                                    binary_comparison= binary_comparison,
                                                    schema_comparison = schema_comparison)
            ARI <- mclust::adjustedRandIndex(truthS, clustering$clustering)

            results <- chen_marginal
            matrix_pvalues <- results$matrix_pvalues %>% data.frame()
            colnames(matrix_pvalues)  <- paste("X", 1:ncol(matrix_pvalues), sep = "")
            truthF_comparison <- results$truthF_comparison %>% data.frame()
            colnames(truthF_comparison) <- paste("X", 1:ncol(truthF_comparison), sep = "")
            comparison <- results$comparison
            time_ = results$time_mean_pvalue

            gene_truehypo <- ((truthF_comparison)) %>%
              data.frame() %>%
              rowid_to_column(var = "variables") %>%
              pivot_longer(cols = c(-variables), names_to = "cluster", values_to = "truehypothesis") %>%
              left_join(
                comparison %>%
                  data.frame() %>%
                  rowid_to_column(var = "cluster") %>%
                  mutate(cluster = paste("X", cluster, sep = "")) %>%
                  mutate(comparison = paste(X1, "vs", X2)) %>%
                  dplyr::select(-c(X1, X2)),
                by = c("cluster")
              ) %>% mutate(variables = as.character(variables))

            pvalues_ext <- matrix_pvalues %>%
              data.frame() %>%
              rownames_to_column(var = "variables") %>%
              pivot_longer(cols = c(-variables),
                           names_to = "cluster",
                           values_to = "pvalues") |>
              left_join(gene_truehypo ,by = c("variables", "cluster"))


            df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                       ndraws = 0, a = a, epsilon = 0,
                                       method_inference = "Chen-marginal",
                                       method_clustering = cl_meth,
                                       exact = TRUE,
                                       sigma = sigma,  rho = rho,
                                       pvalues_ext,
                                       ARI = ARI, time = as.numeric(time_ )))

            rm(chen_marginal)
            rm(gene_truehypo)
            rm(pvalues_ext)
          }
        }
      }


      ################ POCLIN ###################

      if("poclin" %in% methods){
        for(dist in c("rescale")){
          poclin <- pipeline_poclin(X = X, K = K, Sigma = Sigma,
                                    aggregation = dist,
                                    comparison = comparison,
                                    MatMean = MatMean,
                                    binary_comparison= binary_comparison,
                                    schema_comparison = schema_comparison)
          ARI <- mclust::adjustedRandIndex(truthS, poclin$clustering)

          results <- poclin
          matrix_pvalues <- results$matrix_pvalues %>% data.frame()
          colnames(matrix_pvalues)  <- paste("X", 1:ncol(matrix_pvalues), sep = "")
          truthF_comparison <- results$truthF_comparison %>% data.frame()
          colnames(truthF_comparison) <- paste("X", 1:ncol(truthF_comparison), sep = "")
          comparison <- results$comparison
          time_ = results$time_mean_pvalue

          gene_truehypo <- ((truthF_comparison)) %>%
            data.frame() %>%
            rowid_to_column(var = "variables") %>%
            pivot_longer(cols = c(-variables), names_to = "cluster", values_to = "truehypothesis") %>%
            # mutate(nosignal = case_when(variables %in% which(truthF > max(truthS)) ~ "no signal",
            #                             .default = "signal in other clusters")) %>%
            left_join(
              comparison %>%
                data.frame() %>%
                rowid_to_column(var = "cluster") %>%
                mutate(cluster = paste("X", cluster, sep = "")) %>%
                mutate(comparison = paste(X1, "vs", X2)) %>%
                dplyr::select(-c(X1, X2)),
              by = c("cluster")
            ) %>% mutate(variables = as.character(variables))

          pvalues_ext <- matrix_pvalues %>%
            data.frame() %>%
            rownames_to_column(var = "variables") %>%
            pivot_longer(cols = c(-variables),
                         names_to = "cluster",
                         values_to = "pvalues") |>
            left_join(gene_truehypo ,by = c("variables", "cluster"))


          df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                     ndraws = 0, a = a, epsilon = 0,
                                     method_inference = "poclin",
                                     method_clustering = dist,
                                     exact = TRUE,
                                     sigma = sigma,  rho = rho,
                                     pvalues_ext,
                                     ARI = ARI, time = as.numeric(time_ )))
          rm(poclin)
          rm(gene_truehypo)
          rm(pvalues_ext)

        }
      }

      ################ TN TEST ###################
      CLUSTERINGS = c("HAC")
      if("tn.test" %in% methods){
        if (a == 0){
          for(cl_meth in CLUSTERINGS){
            tn_test <- pipeline_tn_test(X = X, K = K,
                                        clustering_method = cl_meth,
                                        comparison = comparison,
                                        MatMean = MatMean,
                                        binary_comparison= binary_comparison,
                                        schema_comparison = schema_comparison)
            ARI <- mclust::adjustedRandIndex(truthS, tn_test$clustering)

            results <- tn_test
            matrix_pvalues <- results$matrix_pvalues %>% data.frame()
            colnames(matrix_pvalues)  <- paste("X", 1:ncol(matrix_pvalues), sep = "")
            truthF_comparison <- results$truthF_comparison %>% data.frame()
            colnames(truthF_comparison) <- paste("X", 1:ncol(truthF_comparison), sep = "")
            comparison <- results$comparison
            time_ = results$time_mean_pvalue

            gene_truehypo <- ((truthF_comparison)) %>%
              data.frame() %>%
              rowid_to_column(var = "variables") %>%
              pivot_longer(cols = c(-variables), names_to = "cluster", values_to = "truehypothesis") %>%
              left_join(
                comparison %>%
                  data.frame() %>%
                  rowid_to_column(var = "cluster") %>%
                  mutate(cluster = paste("X", cluster, sep = "")) %>%
                  mutate(comparison = paste(X1, "vs", X2)) %>%
                  dplyr::select(-c(X1, X2)),
                by = c("cluster")
              ) %>% mutate(variables = as.character(variables))

            pvalues_ext <- matrix_pvalues %>%
              data.frame() %>%
              rownames_to_column(var = "variables") %>%
              pivot_longer(cols = c(-variables),
                           names_to = "cluster",
                           values_to = "pvalues") |>
              left_join(gene_truehypo ,by = c("variables", "cluster"))


            df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                       ndraws = 0, a = a, epsilon = 0,
                                       method_inference = "tn_test",
                                       method_clustering = cl_meth,
                                       exact = TRUE,
                                       sigma = sigma,  rho = rho,
                                       pvalues_ext,
                                       ARI = ARI, time = as.numeric(time_ )))

          }
        }
      }

      return(df)

    })

    df_end <- Reduce(rbind, (res))

    name_file <- sprintf("sim_univariate_set%s_K=%s_a=%s_n=%s_p=%s_nb_experiment=%s_sigma=%s_rho=%s.RDS",
                         setting, K, a, n, p, nb_experiments, sigma, rho)
    file.path <- paste(path, name_file, sep = "/")
    saveRDS(df_end, file = file.path)

  }

}



