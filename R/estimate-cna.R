#' cluster CNAs and resolve CNA state
#' @export
estimateCNA <- function(input_data) {
  
  # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  cna_est <- input_data$tcn
  for (i in seq_len(nrow(input_data$tcn))) {
    # if mutation overlaps CNA region
    if (mean(input_data$tcn[i,]) >= 2) {
        # set CNA to ceiling of max of CNAs
        cna_est[i,] = ceiling(max(input_data$tcn[i,]))
      } else {
        # otherwise, set CNA to floor of min of CNA
        cna_est[i,] = floor(min(input_data$tcn[i,]))
      }
    }
  warning("need to consider if tcn above/below 2 coexist; MCMC to sample cna")
  return(cna_est)
}

testCNA <- function(var=1, min_cluster_thresh=0.11) {
  cna_tmp <- matrix(c(2.1, 2.9, 1.8, 1.9, 1.9, 1.4, 2.2, 3.1, 1.3), nrow = 3)
  # likelihood assuming not copy number alteration
  cna_neutral_likelihood = rowSums(calcCNALikelihood(cna_tmp, 2, var))
  
  icn = ifelse(rowMeans(cna_tmp)>2, ceiling(apply(cna_tmp, 1, max)), floor(apply(cna_tmp, 1, min)))
  icn = rep(icn, ncol(cna_tmp))
  icn = matrix(icn, nrow=nrow(cna_tmp))
  
  cncf_tmp = ifelse((cna_tmp-2)/(icn-2)<=min_cluster_thresh, 0, (cna_tmp-2)/(icn-2))
  cncf_tmp = ifelse(cncf_tmp>1, 1, cncf_tmp)
  
  mu_tmp = (1 - cncf_tmp) * 2 + cncf_tmp * icn
  var_tmp = (1 - cncf_tmp) ^ 2 * var + cncf_tmp ^ 2 * var
  rowSums(calcCNALikelihood(cna_tmp, mu_tmp, var_tmp))
    
  # for (i in seq_len(nrow(cna_tmp))) {
  #   cncf_tmp = 1
  # }
}

calcCNALikelihood <- function(x, mu, var) {
  return(-log(2*pi)/2-log(var^2)/2-(x-mu)^2/(2*var^2))
}

#' estimation of copy number cellular fraction
#' @export
estimateCNCF <- function(input_data, cna_est) {
  cncf_est <- (input_data$tcn - 2) / (cna_est - 2)
  warning("need to consider if tcn above/below 2 coexist; consider MCMC for sampling")
  return(cncf_est)
}


#' assign CNCF to the closest cluster using euclidean distance
reassignCNCF <- function(cncf_est, w_chain) {
  w_mat <- estimateCCFs(w_chain)
  cncf_update <- cncf_est
  for (i in seq_len(nrow(cncf_est))) {
    # print(cncf_est[i,])
    min_dist = 1000000000
    min_row = 1
    for (j in seq_len(nrow(w_mat))) {
      dis = sqrt(sum((cncf_est[i,] - w_mat[j,])^2))
      if (dis < min_dist) {
        min_dist <- dis
        min_row = j
      }
    }
    cncf_update[i,] <- w_mat[min_row,]
  }
  return(cncf_update)
}

#' assign integer copy number using updated CNCF and input tcn
#' max_cn: maximum integer copy number allowed
reassignCNA <- function(cncf_update, tcn, max_cn = 8) {
  cna_update <- cncf_update
  for (i in seq_len(nrow(cncf_update))) {
    min_dist = 100000000000
    icn = 0
    for (j in seq_len(max_cn+1)) {
      j = j - 1
      dis = sqrt(sum((cncf_update[i,] * j + (1-cncf_update[i,]) * 2 - tcn[i,])^2))
      if (dis < min_dist) {
        min_dist <- dis
        icn = j
      }
    }
    cna_update[i,] = icn
  }
  return(cna_update)
}

#' cluster CNAs and resolve CNA state in input_data$y format
#' @export
estimateCNA1 <- function(input_data) {
  
  # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  cna_est <- input_data$y
  
  for (i in seq_len(nrow(input_data$y))) {
    # if mutation overlaps CNA region
    if (any(input_data$overlap[i, ]==1)) {
      cna_index = which(input_data$overlap[i,]==1)
      # if mean CNA across all samples >= 2
      if (mean(input_data$tcn[cna_index,]) >= 2) {
        # set CNA to ceiling of max of CNAs
        cna_est[i,] = ceiling(max(input_data$tcn[cna_index,]))
      } else {
        # otherwise, set CNA to floor of min of CNA
        cna_est[i,] = floor(min(input_data$tcn[cna_index,]))
      }
    } else {
      cna_est[i,] = 2
    }
  }
  warning("need to consider if tcn above/below 2 coexist; MCMC to sample cna")
  return(cna_est)
}


#' estimation of copy number cellular fraction in input_data$y format
#' @export
estimateCNCF1 <- function(input_data, cna_est) {
  cncf_est <- cna_est
  for (i in seq_len(nrow(cna_est))) {
    # if mutation overlaps CNA region
    if (any(input_data$overlap[i, ]==1)) {
      cna_index = which(input_data$overlap[i,]==1)
      cncf_est[i,] = (input_data$tcn[cna_index,] - 2) / (cna_est[i,] - 2)
    } else {
      cncf_est[i,] = 0
    }
  }
  warning("need to consider if tcn above/below 2 coexist; consider MCMC for sampling")
  return(cncf_est)
}