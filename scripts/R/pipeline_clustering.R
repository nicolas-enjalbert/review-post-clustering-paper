pipeline_kmeans <- function(X, K) {
  time_clustering_start <- Sys.time()
  clustering <- kmeans(X, centers = K)$cluster
  time_clustering_stop <- Sys.time()

  # rename clusters as truthS
  clustering <- as.factor(clustering)
  old <- unique((clustering))
  levels(clustering) <- order(old)
  clustering <- as.numeric(as.character(clustering))

  return(list(
    clustering = clustering,
    time_clustering = time_clustering_stop - time_clustering_start
  ))
}

pipeline_HAC <- function(X, K, distance = "euclidean", linkage = "ward.D2") {
  time_clustering_start <- Sys.time()
  clustering <- cutree(hclust(dist(X, method = distance), method = linkage),
    k = K
  )
  time_clustering_stop <- Sys.time()

  # rename clusters as truthS
  clustering <- as.factor(clustering)
  old <- unique((clustering))
  levels(clustering) <- order(old)
  clustering <- as.numeric(as.character(clustering))

  return(list(
    clustering = clustering,
    time_clustering = time_clustering_stop - time_clustering_start
  ))
}

pipeline_GMM <- function(X, K) {
  time_clustering_start <- Sys.time()
  model <- mixmodGaussianModel(listModels = c("Gaussian_p_L_I"))
  res0 <- mixmodCluster(as.data.frame(X), nbCluster = K, models = model)
  clustering <- res0@bestResult@partition
  time_clustering_stop <- Sys.time()

  # rename clusters as truthS
  clustering <- as.factor(clustering)
  old <- unique((clustering))
  levels(clustering) <- order(old)
  clustering <- as.numeric(as.character(clustering))

  return(list(
    clustering = clustering,
    time_clustering = time_clustering_stop - time_clustering_start
  ))
}
