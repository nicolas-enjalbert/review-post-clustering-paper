library(future)
library(future.apply)
library(tidyverse)
library(Rmixmod)
library(Matrix)
library(clusterpval)
library(KmeansInference)
library(mclust)
source("scripts/R/function_yun.R")
source("scripts/R/utils.R")
source("scripts/R/pipeline_clustering.R")
source("scripts/R/pipeline_inference.R")


comparison_multivariate_inference_methods <- function(
    methods = c("naif", "data.thinning", "YFB", "Cond"),
    clusterings = c("HAC", "kmeans", "GMM"),
    nb_experiments = 500, setting = c(1,2),
    ns = 100,
    ps = 3,
    as = c(0,10),
    Ndraws = c(1000),
    epsilon_dts = c(0.1, 0.3, 0.5, 0.7, 0.9),
    sigmas = c(1),
    rhos = c(0),
    type_estimation = c("known", "all", "intra"),
    path =sprintf("results/simulation_multivariate_setting%s", setting),
    nb_workers = min(40, availableCores()/4)
){
  plan(multisession, workers = nb_workers)
  set.seed(130124)

  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  configs <- expand.grid(ns = ns,
                         ps = ps,
                         as = as,
                         sigmas = sigmas,
                         rhos = rhos)
  seq_configs <- 1:nrow(configs)

  if (setting == 1){
    K = 2
  } else if (setting == 2){
    K = 3
  }
  HACMethods <- c("ward.D")

  for(cc in seq_configs){
    config <- configs[cc,]
    a <- config[["as"]]
    n <- config[["ns"]]
    p <- config[["ps"]]
    sigma <- config[["sigmas"]]
    rho <- config[["rhos"]]
    iso = rho == 0

    print(sprintf("a = %s ; n = %s ; p = %s ; sigma = %s ; rho = %s",
                  a, n, p, sigma, rho))


    res <- future.apply::future_lapply(1:nb_experiments, future.seed = TRUE, FUN = function(nexp) {

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


      df <- data.frame()

      ###### HAC
      if("Cond" %in% methods){

        if(rho == 0){
          if("HAC" %in% clusterings){
            for (method in HACMethods){
              cl_HAC <- function(X){
                tree <- hclust(dist(X, method = "euclidean")^2, method = method)
                cutree(tree, k = K)
              }

              clustering_HAC <- cl_HAC(X)
              est_cov <- estimation_cov(X, clustering_HAC, rho)
              for(type_cov in type_estimation){

                sigma_used <- switch(type_cov,
                                     known = sigma,
                                     all = est_cov$sigma_all[1],
                                     intra = est_cov$sigma_intra[1])

                for(ndraws in Ndraws[Ndraws != 0]){
                  t1 <- Sys.time()
                  Gao_MC <- test_clusters_approx(X, k1 = 1, k2 = 2, iso = iso,
                                                 sig = sigma_used, ndraws = ndraws,
                                                 cl_fun = cl_HAC)
                  time_ <- difftime(Sys.time(), t1, units = "secs")
                  ARI <- mclust::adjustedRandIndex(Gao_MC$clusters, truthS)
                  df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                             ndraws = ndraws, a = a, epsilon = 0,
                                             method_inference = "Gao",
                                             method_clustering = paste("HAC", method,
                                                                       sep = "_"),
                                             exact = FALSE,
                                             sigma = sigma,  rho = rho,
                                             type_estimation = type_cov,
                                             pvaleur = Gao_MC$pval,
                                             stat = Gao_MC$stat,
                                             ARI = ARI, time = as.numeric(time_ )))
                }

                if(0 %in% Ndraws){
                  t1 <- Sys.time()
                  tree <- hclust(dist(X, method = "euclidean")^2, method = method)
                  Gao_exact <- test_hier_clusters_exact(X = X, link = method,
                                                        hcl = tree, K = K, k1 = 1,
                                                        k2 = 2, iso = iso, sig = sigma_used)
                  time_ <- difftime(Sys.time(), t1, units = "secs")
                  ARI <- mclust::adjustedRandIndex(cutree(tree, k = K), truthS)
                  df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                             a = a, epsilon = 0,
                                             method_inference = "Gao",
                                             method_clustering = paste("HAC", method, sep = "_"),
                                             exact = TRUE,
                                             sigma = sigma,  rho = rho,
                                             type_estimation = type_cov,
                                             pvaleur = Gao_exact$pval,
                                             stat = Gao_exact$stat,
                                             ARI = ARI, time = as.numeric(time_ )))

                }
              }
            }
          }

          ######## kmeans
          if("kmeans" %in% clusterings){
            seed <- sample(1:10000, 1)
            cl_kmeans <- function(X){
              set.seed(seed)
              kmeans(x = X, centers = K, algorithm = "Lloyd", iter.max = 1000)$cluster
            }

            clustering_KM <- cl_kmeans(X)
            est_cov <- estimation_cov(X, clustering_KM, rho)

            for(type_cov in type_estimation){

              sigma_used <- switch(type_cov,
                                   known = sigma,
                                   all = est_cov$sigma_all[1],
                                   intra = est_cov$sigma_intra[1])

              for(ndraws in Ndraws[which(Ndraws!=0)]){
                t1 <- Sys.time()
                kmeans_MC <- test_clusters_approx(X, k1 = 1, k2 = 2, iso = iso, sig = sigma_used,
                                                  ndraws = ndraws, cl_fun = cl_kmeans)
                time_ <- difftime(Sys.time(), t1, units = "secs")
                ARI <- mclust::adjustedRandIndex(kmeans_MC$clusters, truthS)
                df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = ndraws,
                                           a = a, epsilon = 0,
                                           method_inference = "Gao",
                                           method_clustering = "kmeans",
                                           exact = FALSE,
                                           sigma = sigma,  rho = rho,
                                           type_estimation = type_cov,
                                           pvaleur = kmeans_MC$pval,
                                           stat = kmeans_MC$stat,
                                           ARI = ARI, time = as.numeric(time_ )))
              }

              if(0 %in% Ndraws){
                t1 <- Sys.time()
                kmeans_exact <- kmeans_inference(X = X, k = K, cluster_1 = 1, cluster_2 = 2,
                                                 iso = iso, sig = sigma_used, iter.max = 1000,
                                                 seed = seed, verbose = TRUE)
                time_ <- difftime(Sys.time(), t1, units = "secs")
                ARI <- mclust::adjustedRandIndex(kmeans_exact$final_cluster, truthS)
                df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                           a = a, epsilon = 0,
                                           method_inference = "Gao",
                                           method_clustering = "kmeans",
                                           exact = TRUE,
                                           sigma = sigma,  rho = rho,
                                           type_estimation = type_cov,
                                           pvaleur = kmeans_exact$pval,
                                           stat = kmeans_exact$test_stat,
                                           ARI = ARI, time = as.numeric(time_ )))
              }
            }
          }

          ####### GMM
          if("GMM" %in% clusterings){
            GMMclustering <- function(X){
              model <- mixmodGaussianModel(listModels = c("Gaussian_p_L_I"))
              res0 <- mixmodCluster(as.data.frame(X), nbCluster = K, models = model)
              clustering <- res0@bestResult@partition
              return(clustering)
            }

            clustering_GMM <- GMMclustering(X)
            est_cov <- estimation_cov(X, clustering_GMM, rho)

            for(type_cov in type_estimation){

              sigma_used <- switch(type_cov,
                                   known = sigma,
                                   all = est_cov$sigma_all[1],
                                   intra = est_cov$sigma_intra[1])

              for(ndraws in Ndraws[Ndraws != 0]){
                t1 <- Sys.time()
                GMM_MC = test_clusters_approx(X, k1 = 1, k2 = 2, iso = iso,
                                              sig = sigma_used,
                                              ndraws = ndraws,
                                              cl_fun = GMMclustering)
                time_ <- difftime(Sys.time(), t1, units = "secs")
                ARI <- mclust::adjustedRandIndex(GMM_MC$clusters, truthS)
                df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p,
                                           ndraws = ndraws,
                                           a = a, epsilon = 0,
                                           method_inference = "Gao",
                                           method_clustering = "GMM",
                                           exact = FALSE,
                                           sigma = sigma,  rho = rho,
                                           type_estimation = type_cov,
                                           pvaleur = GMM_MC$pval,
                                           stat = GMM_MC$stat,
                                           ARI = ARI, time = as.numeric(time_ )))
              }
            }
          }
        } else {
        if("HAC" %in% clusterings){
          cl_HAC <- function(X){
            tree <- hclust(dist(X, method = "euclidean")^2, method = "ward.D")
            cutree(tree, k = K)
          }

          clustering_HAC <- cl_HAC(X)
          est_cov <- estimation_cov(X, clustering_HAC, rho)
          for(type_cov in type_estimation){

            Sigma_used <- switch(type_cov,
                                 known = Sigma,
                                 all = est_cov$Sigma_all,
                                 intra = est_cov$Sigma_intra)


            U <- matrixNormal::I(n) # U: dependence between observations
            t1 <- Sys.time()
            javi_exact = test.clusters.hc(X, U, Sigma_used, NC = 2, clusters = c(1,2),
                                          plot = FALSE, linkage = "ward.D")
            time_ <- difftime(Sys.time(), t1, units = "secs")
            ARI <- mclust::adjustedRandIndex(javi_exact$hcl, truthS)
            df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                       a = a, epsilon = 0,
                                       method_inference = "general_dependance",
                                       method_clustering = "HAC",
                                       exact = TRUE,
                                       sigma = sigma,  rho = rho,
                                       type_estimation = type_cov,
                                       pvaleur = javi_exact$pvalue,
                                       stat = javi_exact$stat,
                                       ARI = ARI, time = as.numeric(time_ )))
          }
        }
      }}


      ####### yun HAC ward

      if("YFB" %in% methods){
        if(0 %in% Ndraws){
        t1 <- Sys.time()
        yun_exact = fun_proposed_exact(X, cl_meth = "ward")
        time_ <- difftime(Sys.time(), t1, units = "secs")
        ARI <- mclust::adjustedRandIndex(yun_exact$cl, truthS)
        df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                   a = a, epsilon = 0,
                                   method_inference = "unknown_variance",
                                   method_clustering = "HAC",
                                   exact = TRUE,
                                   sigma = sigma,  rho = rho,
                                   type_estimation = "",
                                   pvaleur = yun_exact$pval,
                                   stat = 0,
                                   ARI = ARI, time = as.numeric(time_ )))
        }

        alpha = 0.05
        for(ndraws in Ndraws[Ndraws != 0]){
          t1 <- Sys.time()
          yun_ISMC = fun_proposed_approx(X, K = K, k1 = 1, k2 = 2, ndraws = ndraws,
                                         alpha = alpha, cl_meth = "ward")
          time_ <- difftime(Sys.time(), t1, units = "secs")
          ARI <- mclust::adjustedRandIndex(yun_ISMC$cl, truthS)
          df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = ndraws,
                                     a = a, epsilon = 0,
                                     method_inference = "unknown_variance",
                                     method_clustering = "HAC",
                                     exact = FALSE,
                                     sigma = sigma,  rho = rho,
                                     type_estimation = "",
                                     pvaleur = yun_ISMC$pval,
                                     stat = 0,
                                     ARI = ARI, time = as.numeric(time_ )))

        }
      }



      ####### naif

      if("naif" %in% methods){
        if(rho == 0 & "known" %in% type_estimation){

          if("kmeans" %in% clusterings){
            t1 <- Sys.time()
            clustering <- kmeans(X, centers = K)$cluster
            eta <- eta_comparison_mean(clustering = clustering, k1 = 1, k2 = 2)
            naif_kmeans <- multivariate_Z_test(eta= eta, X = X, sig = sigma)
            time_ <- difftime(Sys.time(), t1, units = "secs")
            ARI <- mclust::adjustedRandIndex(clustering, truthS)
            df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                       a = a, epsilon = 0,
                                       method_inference = "naif",
                                       method_clustering = "kmeans",
                                       exact = TRUE,
                                       sigma = sigma,  rho = rho,
                                       type_estimation = "known",
                                       pvaleur = naif_kmeans$pval,
                                       stat = naif_kmeans$stat,
                                       ARI = ARI, time = as.numeric(time_ )))
          }

          if("HAC" %in% clusterings){
            t1 <- Sys.time()
            clustering <- cutree(hclust(dist(X), method = "ward.D2"), k = K)
            eta <- eta_comparison_mean(clustering = clustering, k1 = 1, k2 = 2)
            naif_HAC <- multivariate_Z_test(eta= eta, X = X, sig = sigma)
            time_ <- difftime(Sys.time(), t1, units = "secs")
            ARI <- mclust::adjustedRandIndex(clustering, truthS)
            df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                       a = a, epsilon = 0,
                                       method_inference = "naif",
                                       method_clustering = "HAC",
                                       exact = TRUE,
                                       sigma = sigma,  rho = rho,
                                       type_estimation = "known",
                                       pvaleur = naif_HAC$pval,
                                       stat = naif_HAC$stat,
                                       ARI = ARI, time = as.numeric(time_ )))
          }
          if("GMM" %in% clusterings){
            t1 <- Sys.time()
            model <- mixmodGaussianModel(listModels = c("Gaussian_p_L_I"))
            res0 <- mixmodCluster(as.data.frame(X), nbCluster = K, models = model)
            clustering <- res0@bestResult@partition
            eta <- eta_comparison_mean(clustering = clustering, k1 = 1, k2 = 2)
            naif_GMM <- multivariate_Z_test(eta= eta, X = X, sig = sigma)
            time_ <- difftime(Sys.time(), t1, units = "secs")
            ARI <- mclust::adjustedRandIndex(clustering, truthS)
            df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                       a = a, epsilon = 0,
                                       method_inference = "naif",
                                       method_clustering = "GMM",
                                       exact = TRUE,
                                       sigma = sigma,  rho = rho,
                                       type_estimation = "known",
                                       pvaleur = naif_GMM$pval,
                                       stat = naif_GMM$stat,
                                       ARI = ARI, time = as.numeric(time_ )))
          }
        } else {
          print("Naive method can be only used with rho = 0 and Sigma known")
        }
      }


      ######### Data thinning

      if("data.thinning" %in% methods){
        if("known" %in% type_estimation){
          for (epsilon_dt in epsilon_dts){

            t1_dt <- Sys.time()
            XDT <- datathin_gaussian_multivariate(X, epsilon = epsilon_dt, Sigma = Sigma)
            Xclustering <- XDT$Xtr
            Xtest <- XDT$Xte
            time_dt <- difftime(Sys.time(), t1_dt, units = "secs")

            if("kmeans" %in% clusterings){
              t1 <- Sys.time()
              clustering <- kmeans(Xclustering, centers = K)$cluster
              eta <- eta_comparison_mean(clustering = clustering, k1 = 1, k2 = 2)
              naif_kmeans <- multivariate_Z_test(eta= eta, X = Xtest, sig = sqrt(1-epsilon_dt)*sigma)
              time_ <- difftime(Sys.time(), t1, units = "secs")
              ARI <- mclust::adjustedRandIndex(clustering, truthS)
              df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                         a = a, epsilon = epsilon_dt,
                                         method_inference = "datathinning",
                                         method_clustering = "kmeans",
                                         exact = TRUE,
                                         sigma = sigma,  rho = rho,
                                         type_estimation = "known",
                                         pvaleur = naif_kmeans$pval,
                                         stat = naif_kmeans$stat,
                                         ARI = ARI, time = as.numeric(time_ + time_dt )))
            }

            if("HAC" %in% clusterings){
              t1 <- Sys.time()
              clustering <- cutree(hclust(dist(Xclustering), method = "ward.D2"), k = K)
              eta <- eta_comparison_mean(clustering = clustering, k1 = 1, k2 = 2)
              naif_HAC <- multivariate_Z_test(eta= eta, X = Xtest, sig = sqrt(1-epsilon_dt)*sigma)
              time_ <- difftime(Sys.time(), t1, units = "secs")
              ARI <- mclust::adjustedRandIndex(clustering, truthS)
              df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                         a = a, epsilon = epsilon_dt,
                                         method_inference = "datathinning",
                                         method_clustering = "HAC",
                                         exact = TRUE,
                                         sigma = sigma,  rho = rho,
                                         type_estimation = "known",
                                         pvaleur = naif_HAC$pval,
                                         stat = naif_HAC$stat,
                                         ARI = ARI, time = as.numeric(time_ + time_dt )))
            }

            if("GMM" %in% clusterings){
              t1 <- Sys.time()
              model <- mixmodGaussianModel(listModels = c("Gaussian_p_L_I"))
              res0 <- mixmodCluster(as.data.frame(Xclustering), nbCluster = K, models = model)
              clustering <- res0@bestResult@partition
              eta <- eta_comparison_mean(clustering = clustering, k1 = 1, k2 = 2)
              naif_GMM <- multivariate_Z_test(eta= eta, X = Xtest, sig = sqrt(1-epsilon_dt)*sigma)
              time_ <- difftime(Sys.time(), t1, units = "secs")
              ARI <- mclust::adjustedRandIndex(clustering, truthS)
              df <- rbind(df, data.frame(nexp = nexp, K = K, n = n, p = p, ndraws = 0,
                                         a = a, epsilon = epsilon_dt,
                                         method_inference = "datathinning",
                                         method_clustering = "GMM",
                                         exact = TRUE,
                                         sigma = sigma,  rho = rho,
                                         type_estimation = "known",
                                         pvaleur = naif_GMM$pval,
                                         stat = naif_GMM$stat,
                                         ARI = ARI, time = as.numeric(time_ + time_dt)))
            }

          }
        } else {
          print("Data thinning method only works with known covariance matrix.")
        }
      }
      return(df)
    }
    )
    df_end <- Reduce(rbind, (res))

    name_file <- sprintf("sim_multi_set%s_K=%s_a=%s_n=%s_p=%s_nb_experiment=%s_sigma=%s_rho=%s.RDS",
                         setting, K, a, n, p, nb_experiments, sigma, rho)
    file.path <- paste(path, name_file, sep = "/")
    saveRDS(df_end, file = file.path)
  }




}
