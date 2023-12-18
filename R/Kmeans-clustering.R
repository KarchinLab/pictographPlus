#' run kmeans clusterig to cluster mutations/CNAs with similar mcf/cncf
#' @export
KmeansClustering <- function(total_est, total_sep_list, Kmax=5, cluster_mutation_thresh=1, cluster_diff_thresh=0.05) {
  clusters <- list()
  for (index in seq_len(length(total_sep_list))) {
    df = as.data.frame(total_est[total_sep_list[[index]],])
    sample_index = which(strsplit(names(total_sep_list[index]), split="")[[1]]=="1")
    df1 = df[, sample_index, drop=F]
    bestK = 1
    bestCluster = c()
    toBreak = F
    # loop through all possible k until hits Kmax, or error in kmeans, or encounter cluster with one one mutation, 
    # or encounter two clusters where the mean mcf/cncf is less than certain threshold in all samples conditioned on sample presence (ignore columns with all zeros)
    for (k in seq_len(Kmax)) {
      tryCatch({
        km.res = kmeans(df1, k)
      }, error = function(e) {
        toBreak <<- T
      })
      if (toBreak) {
        break
      }
      
      if (k>1) {
        # check whether all clusters contain more than cluster_mutation_thresh mutations
        for (i in seq_len(k)) {
          if (length(which(km.res$cluster==i))<=cluster_mutation_thresh) {
            toBreak = T
          }
        }
        # check whether mean difference between any two clusters is less than cluster_diff_thresh among all samples
        if (!toBreak) {
          cluster_mean <- list()
          for (i in seq_len(k-1)) {
            for (j in seq(i+1, k)) {
              diff = abs(colMeans(df1[which(km.res$cluster==i),]) - colMeans(df1[which(km.res$cluster==j),]))
              if (all(diff<=cluster_diff_thresh)) {
                toBreak = T
              }
            }
          }
        }
      }
      if (toBreak) {
        break
      }
      bestK = k
      bestCluster = km.res$cluster
    }
    clusters[[names(total_sep_list[index])]] = bestCluster
  }
  return(clusters)
}