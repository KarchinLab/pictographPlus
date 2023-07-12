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


#' estimation of copy number cellular fraction
#' @export
estimateCNCF <- function(input_data, cna_est) {
  cncf_est <- (input_data$tcn - 2) / (cna_est - 2)
  warning("need to consider if tcn above/below 2 coexist; consider MCMC for sampling")
  return(cncf_est)
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