#' estimate cluster mcf/cncf by taking the mean
#' @export
estimateClusterMCF <- function(total_est, clusters) {
  cluster = 1
  clusterMCF = matrix(,ncol=ncol(total_est))
  colnames(clusterMCF) = colnames(total_est)
  rowname = c()
  for (index in seq_len(length(clusters))) {
    for (k in unique(clusters[[index]])) {
      clusterMCF_tmp = total_est[names(which(clusters[[index]]==k)),]
      if (nrow(clusterMCF_tmp)>1) {
        clusterMCF_tmp = colMeans(clusterMCF_tmp)
      }
      clusterMCF = rbind(clusterMCF, clusterMCF_tmp)
      rowname <- append(rowname, cluster)
      cluster = cluster + 1
    }
  }
  clusterMCF <- clusterMCF[-c(1),,drop=F]
  rownames(clusterMCF) = rowname
  return(clusterMCF)
}

#' get cluster assignment
#' @export
getClusterAssignment <- function(total_est, clusters) {
  cluster = 1
  clusterAssignment = matrix(,ncol=2)
  colnames(clusterAssignment) = c('Mut_ID', 'Cluster')
  for (index in seq_len(length(clusters))) {
    for (k in unique(clusters[[index]])) {
      clusterAssignment_tmp = matrix(,nrow=length(names(which(clusters[[index]]==k))), ncol=2)
      clusterAssignment_tmp[,1] = names(which(clusters[[index]]==k))
      clusterAssignment_tmp[,2] = cluster
      clusterAssignment = rbind(clusterAssignment, clusterAssignment_tmp)
      cluster = cluster + 1
    }
  }
  clusterAssignment <- noquote(clusterAssignment[-c(1),,drop=F])
  return(clusterAssignment)
}
