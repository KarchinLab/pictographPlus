#' Determine the most probable cluster CCF values by taking the mode of the posterior distributions
#' 
#' @export
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @return matrix of estimated cluster CCFs
estimateCCFs <- function(w_chain) {
  S <- numberSamples(w_chain)
  K <- numberClusters(w_chain)
  # density plot 
  w.dens <- ggplot(w_chain, aes(x = value)) +
    geom_density() +
    facet_wrap(~Parameter, ncol = S, scales = "free_y") +
    theme_light()
  # find peak for MAP w
  w.dens.p <- ggplot_build(w.dens)$data[[1]]
  w.map <- w.dens.p %>%
    as_tibble() %>%
    group_by(PANEL) %>%
    summarize(value = x[max(y) == y])
  w.map <- w.map %>%
    mutate(Parameter = unique(w_chain$Parameter),
           value_rounded = round(value, 2))
  # return w matrix
  w.map.matrix <- matrix(w.map$value_rounded, K, S, byrow=TRUE)
  return(w.map.matrix)
}

#' @importFrom stringr str_replace
numberSamples <- function(mcf_stats){
  params <- as.character(mcf_stats$Parameter)    
  nSamples <- strsplit(params, ",") %>%
    sapply("[", 2) %>%
    stringr::str_replace("\\]", "") %>%
    as.numeric() %>%
    max()
  nSamples
}

#' @importFrom stringr str_replace
numberClusters <- function(mcf_stats){
  params <- as.character(mcf_stats$Parameter)
  K <- strsplit(params, ",") %>%
    sapply("[", 1) %>%
    str_replace("w\\[", "") %>%
    as.numeric() %>%
    max()
  K
}

#' Determine the most probable cluster CCF values by taking the mode of the posterior distributions
#' 
#' @export
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @param Sample_ID Vector of sample IDs, same order as provided as input data (e.g. indata$Sample_ID)
#' @return A tibble of estimated cluster CCFs in each sample 
writeClusterCCFsTable <- function(w_chain, Sample_ID = NULL) {
  map_w <- as.data.frame(estimateCCFs(w_chain))
  
  if (is.null(Sample_ID)) {
    Sample_ID <- paste0("Sample ", 1:ncol(map_w))
  }
  colnames(map_w) <- Sample_ID
  map_w <- map_w %>%
    as_tibble() %>%
    bind_cols(tibble(Cluster = 1:nrow(map_w)), .)
  return(map_w)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability
#' 
#' @export
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
estimateClusterAssignments <- function(z_chain) {
  it <- max(z_chain$Iteration)
  mcmc_z <- z_chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              maxiter=it) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(value=value[probability==max(probability)])
  
  # choose first cluster if equal probability
  map_z_count <- map_z %>% 
    group_by(Parameter) %>%
    summarize(map_count = n()) %>%
    ungroup()
  if (any(map_z_count$map_count > 1)) {
    mut_ind <- which(map_z_count$map_count > 1)
    for (i in mut_ind) {
      dup_var <- as.numeric(gsub("z\\[|]", "", map_z_count$Parameter[i]))
      map_z_dups <- which(gsub("z\\[|]", "", map_z$Parameter) == dup_var)
      dup_ind <- map_z_dups[-1]
      map_z <- map_z[-dup_ind, ]
    }
  }
  return(map_z)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability. 
#' 
#' @export
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
#' @param Mut_ID Vector of mutation IDs, same order as provided as input data (e.g. indata$Mut_ID)
#' @return A tibble listing mutation IDs and their cluster assignments
writeClusterAssignmentsTable <- function(z_chain, cncf=NULL, Mut_ID = NULL) {
  map_z <- estimateClusterAssignments(z_chain) 
  if (is.null(Mut_ID)) {
    Mut_ID <- paste0("Mut", 1:nrow(map_z))
  }
  map_z <- map_z %>%
    mutate(Mut_ID = Mut_ID, Cluster = value) %>%
    select(Mut_ID, Cluster) %>%
    arrange(Cluster)
  
  if (!is.null(cncf)) {
    for (i in seq_len(nrow(cncf))) {
      cls = which(apply(w_mat, 1, function(x) return(all(x == cncf_update[i,]))))
      map_z <- map_z %>% add_row(Mut_ID=rownames(cncf)[i], Cluster=cls)
    }
  }
  
  map_z <- map_z %>%
    arrange(Cluster)
  return(map_z)
}

#' Collect chains for best K of each mutation set 
#' 
#' @export
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param chosen_K (Optional) Vector of K to choose for each mutation set, in the same order as all_set_results. If left blank, function will select best K automatically selected by \code{clusterSep}
collectBestKChains <- function(all_set_results, chosen_K = NULL) {
  # best_set_chains <- lapply(all_set_results, function(x) x$all_chains[[length(x$all_chains)]])
  if (is.null(chosen_K)) {
    best_set_chains <- lapply(all_set_results, function(x) x$best_chains)
  } else {
    best_set_chains <- mapply(function(set_res, choose_K) set_res$all_chains[[choose_K]],
                              set_res = all_set_results,
                              chosen_K,
                              SIMPLIFY = FALSE)
  }
  return(best_set_chains)
}

#' Relabel chains for all sets and merge 
#' 
#' @export
#' @import dplyr
#' @param best_set_chains List of lists of MCMC chains (w_chain, z_chain, ystar_chain) for each mutation set
#' @param indata List of input data objects (same as provided to clusterSep)
mergeSetChains <- function(best_set_chains, indata) {
  best_K_vals <- unname(sapply(best_set_chains, function(x) max(x$z_chain$value)))
  sep_list <- separateMutationsBySamplePresence(indata)
  
  # first set doesn't need to change cluster labels
  w_chain <- best_set_chains[[1]]$w_chain
  temp_z_chain <- best_set_chains[[1]]$z_chain
  temp_ystar_chain <- best_set_chains[[1]]$ystar_chain
  
  if (length(best_set_chains) > 1) {
    # still need to change mutation indices if more than 1 box
    z_chain <- relabel_z_chain_mut_only(temp_z_chain, sep_list[[1]]$mutation_indices)
    ystar_chain <- relabel_ystar_chain(temp_ystar_chain,
                                       sep_list[[1]]$mutation_indices)
    for (i in 2:length(best_set_chains)) {
      temp_w_chain <- best_set_chains[[i]]$w_chain
      temp_z_chain <- best_set_chains[[i]]$z_chain
      temp_ystar_chain <- best_set_chains[[i]]$ystar_chain
      new_cluster_labels <- seq_len(best_K_vals[i]) + sum(best_K_vals[1:(i-1)])
      
      temp_relabeled_w_chain <- relabel_w_chain(temp_w_chain, new_cluster_labels)
      temp_relabeled_z_chain <- relabel_z_chain(temp_z_chain, new_cluster_labels, 
                                                sep_list[[i]]$mutation_indices)
      temp_relabeled_ystar_chain <- relabel_ystar_chain(temp_ystar_chain,
                                                        sep_list[[i]]$mutation_indices)
      
      w_chain <- rbind(w_chain, temp_relabeled_w_chain)
      z_chain <- rbind(z_chain, temp_relabeled_z_chain)
      ystar_chain <- rbind(ystar_chain, temp_relabeled_ystar_chain)
    }
  } else {
    z_chain <- temp_z_chain
    ystar_chain <- temp_ystar_chain
  }
  
  # set levels for Parameter
  w_chain <- w_chain %>% 
    mutate(k = as.numeric(gsub("w\\[", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][2])))) %>%
    arrange(k, s) %>%
    mutate(Parameter = factor(Parameter, levels = unique(w_chain$Parameter))) %>%
    select(Iteration, Chain, Parameter, value)
  
  z_chain_param_order <- tibble(Parameter = unique(z_chain$Parameter)) %>%
    mutate(Variant = as.numeric(gsub("z\\[", "", 
                                     gsub("\\]", "", 
                                          unique(z_chain$Parameter))))) %>%
    arrange(Variant)
  z_chain <- z_chain %>%
    mutate(Parameter = factor(Parameter, levels = z_chain_param_order$Parameter))
  
  ystar_chain <- ystar_chain %>%
    mutate(Mutation_index = as.numeric(gsub("ystar\\[", "",
                                            sapply(ystar_chain$Parameter,
                                                   function(x) strsplit(as.character(x), ",")[[1]][1]))),
           s = as.numeric(gsub("\\]", "",
                               sapply(ystar_chain$Parameter,
                                      function(x) strsplit(as.character(x), ",")[[1]][2]))))
  ystar_chain <- ystar_chain %>%
    arrange(Mutation_index, s) %>%
    mutate(Parameter = factor(Parameter, levels = unique(Parameter)))
  
  chains <- list(w_chain = w_chain,
                 z_chain = z_chain,
                 ystar_chain = ystar_chain)
  return(chains)
}

relabel_z_chain <- function(z_chain, new_cluster_labels, mutation_indices) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  if (length(mutation_indices) != length(unique(z_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in z_chain")
  }
  ## would break when no mutation is assigned to a cluster
  ## poor choice of k, would prob lower the k
  if (length(new_cluster_labels) < length(unique(z_chain$value))) {
    stop("number of supplied new cluster labels does not match the number of clusters in z_chain")
  }
  new_z <- z_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("z\\[", "", 
                                    sapply(z_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_z <- new_z %>%
    mutate(new_i = mutation_indices[i],
           value = new_cluster_labels[new_z$value]) %>%
    mutate(Parameter = paste0("z[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_z)
}

relabel_z_chain_mut_only <- function(z_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  # cluster labels are left unchanged 
  if (length(mutation_indices) != length(unique(z_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in z_chain")
  }
  new_z <- z_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("z\\[", "", 
                                    sapply(z_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_z <- new_z %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("z[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_z)
}

relabel_ystar_chain <- function(ystar_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  i_s <- gsub("ystar\\[|]", "", ystar_chain$Parameter)
  i <- sapply(i_s, function(x) strsplit(x, ",")[[1]][1]) %>%
    as.numeric
  s <- sapply(i_s, function(x) strsplit(x, ",")[[1]][2]) %>%
    as.numeric
  new_ystar <- ystar_chain %>%
    mutate(i = i,
           s = s)
  new_ystar <- new_ystar %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("ystar[", new_i, ",", s, "]")) %>%
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_ystar)
}

relabel_w_chain <- function(w_chain, new_cluster_labels) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  new_w <- w_chain %>% 
    mutate(k = as.numeric(gsub("w\\[", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][2]))))
  if (length(new_cluster_labels) != length(unique(new_w$k))) {
    stop("number of supplied new cluster labels does not match the number of clusters in w_chain")
  }
  new_w <- new_w %>% 
    mutate(k_new = new_cluster_labels[new_w$k]) %>%
    mutate(Parameter = paste0("w[", k_new, ",", s, "]")) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_w)
}
