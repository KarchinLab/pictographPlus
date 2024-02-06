#' run MCMC for all mutations by sample presence
#' @export
runMCMCForAllBoxes <- function(sep_list, max_K = 5, min_mutation_per_cluster = 2, cluster_diff_thresh=0.05,
                               n.iter = 5000, n.burn = 1000, thin = 10,
                               mc.cores = 4, model_type = "spike_and_slab", beta.prior = FALSE, drop_zero = TRUE,
                               inits = list(".RNG.name" = "base::Wichmann-Hill",
                                            ".RNG.seed" = 123)){
  all_set_results <- vector("list", length(sep_list))
  names(all_set_results) <- names(sep_list)
  params = c("z", "w", "ystar")

  for (i in seq_len(length(sep_list))) {
    # i = 2
    temp_box <- sep_list[[i]]
    # Max number of clusters cannot be more than number of mutations/min_mutation_per_cluster
    temp_max_K <- min(max_K, floor(length(temp_box$mutation_indices)/min_mutation_per_cluster))
    temp_max_K <- max(temp_max_K, 1)
    temp_samps_list <- runMutSetMCMC(temp_box, 
                                     n.iter = n.iter, n.burn = n.burn, thin = thin, 
                                     mc.cores = mc.cores,
                                     inits = inits,
                                     temp_max_K = temp_max_K,
                                     model_type = model_type,
                                     params = params,
                                     beta.prior = beta.prior,
                                     drop_zero = drop_zero,
                                     min_mutation_per_cluster = min_mutation_per_cluster, 
                                     cluster_diff_thresh = cluster_diff_thresh)
    
    all_set_results[[i]] <- temp_samps_list
    # break
  }
  return(all_set_results)
}

runMCMC <- function(box_input_data, K, jags.file, inits, params,
                    n.iter=10000, thin=10, n.chains=1,
                    n.adapt=1000, n.burn=1000,
                    beta.prior=FALSE) {
  # if (K > 1) data$K <- K
  if (K > 1) box_input_data$K <- K
  jags.m <- jags.model(jags.file,
                       box_input_data,
                       n.chains = n.chains,
                       inits = inits,
                       n.adapt = n.adapt)
  if (n.burn > 0) update(jags.m, n.burn)
  samps <- coda.samples(jags.m, params, n.iter=n.iter, thin=thin)
  
  if (beta.prior & K > 1) {
    # use initial MCMC to estimate beta priors for identified clusters
    initial_chains <- formatChains(samps)
    est_ccfs <- initial_chains$w_chain %>%
      estimateCCFs %>%
      as.data.frame %>%
      magrittr::set_colnames(1:ncol(.)) %>%
      tibble %>%
      mutate(cluster = 1:nrow(.)) %>%
      pivot_longer(cols = colnames(.)[colnames(.) != "cluster"], 
                   names_to = "sample", 
                   values_to = "ccf") %>%
      mutate(sample = as.numeric(sample))
    
    
    beta_shapes <- lapply(est_ccfs$ccf, estimateBetaPriors)
    cluster_beta_params <- est_ccfs %>%
      mutate(shape1 = sapply(beta_shapes, function(x) x[1]),
             shape2 = sapply(beta_shapes, function(x) x[2]))
    
    cluster_shape1 <- cluster_beta_params %>%
      select(cluster, shape1, sample) %>%
      pivot_wider(names_from = sample,
                  values_from = shape1) %>%
      select(-c(cluster)) %>%
      as.matrix
    cluster_shape2 <- cluster_beta_params %>%
      select(cluster, shape2, sample) %>%
      pivot_wider(names_from = sample,
                  values_from = shape2) %>%
      select(-c(cluster)) %>%
      as.matrix
    
    # run second MCMC with specified beta priors
    box_input_data$cluster_shape1 <- cluster_shape1
    box_input_data$cluster_shape2 <- cluster_shape2
    
    extdir <- system.file("extdata", package="pictograph")
    jags.file.beta <- file.path(extdir, "model-simple-set-beta.jags")
    
    jags.m.beta <- rjags::jags.model(jags.file.beta,
                                     box_input_data,
                                     n.chains = n.chains,
                                     inits = inits,
                                     n.adapt = n.adapt)
    if (n.burn > 0) update(jags.m.beta, n.burn)
    samps.beta <- rjags::coda.samples(jags.m.beta, params, n.iter=n.iter, thin=thin)
    return(samps.beta)
  }
  
  return(samps)
}

reverseDrop <- function(samps, pattern, n.iter) {
  total_sample = nchar(pattern)
  sample_list = vector()
  for (j in seq_len(nchar(pattern))) {
    if (strsplit(pattern, "")[[1]][j] == "1") {
      sample_list <- append(sample_list, j)
    }
  }
  k_list = vector()
  ystar_list = vector()
  # replace current sample id by true sample id from pattern
  for (i in seq_len(length(colnames(samps[[1]])))) {
    # print(colnames(samps[[1]])[i])
    if (startsWith(colnames(samps[[1]])[i], "w")) {
      para <- str_extract_all(colnames(samps[[1]])[i], "[0-9]+")[[1]]
      colnames(samps[[1]])[i] <- paste("w[", para[1], ",", sample_list[strtoi(para[2])], "]", sep = "")
      k_list <- c(k_list, para[1])
    } else if (startsWith(colnames(samps[[1]])[i], "ystar")) {
      para <- str_extract_all(colnames(samps[[1]])[i], "[0-9]+")[[1]]
      colnames(samps[[1]])[i] <- paste("ystar[", para[1], ",", sample_list[strtoi(para[2])], "]", sep = "")
      ystar_list <- c(ystar_list, para[1])
    }
  }
  k_list <- unique(k_list)
  ystar_list <- unique(ystar_list)
  
  # add back dropped samples
  absent_sample <- vector()
  for (sample in seq_len(total_sample)) {
    if (! sample %in% sample_list) {
      absent_sample <- append(absent_sample, sample)
    }
  }
  for (k in seq_len(length(k_list))) {
    for (j in seq_len(length(absent_sample))) {
      col = paste("w[", k_list[k], ",", absent_sample[j], "]", sep = "")
      samps[[1]] <- cbind(samps[[1]], col=0)
      colnames(samps[[1]])[colnames(samps[[1]]) == 'col'] <- col
    }
  }
  for (ystar in seq_len(length(ystar_list))) {
    for (j in seq_len(length(absent_sample))) {
      col = paste("ystar[", ystar_list[ystar], ",", absent_sample[j], "]", sep = "")
      samps[[1]] <- cbind(samps[[1]], col=0)
      colnames(samps[[1]])[colnames(samps[[1]]) == 'col'] <- col
    }
  }
  
  samps[[1]] <- samps[[1]][,order(colnames(samps[[1]]))]
  
  return(samps)
}

getBoxInputData <- function(box) {
  sample_list = vector()
  for (j in 1:ncol(box$y)) {
    if (strsplit(box$pattern, "")[[1]][j] == "1") {
      sample_list <- append(sample_list, j)
    }
  }
  # print(temp_box)
  box_input_data <- list(I = nrow(box$y),
                         S = length(sample_list),
                         y = box$y[,sample_list,drop=FALSE],
                         n = box$n[,sample_list,drop=FALSE],
                         m = box$m[,sample_list,drop=FALSE],
                         icn = box$tcn[,sample_list,drop=FALSE], # integer copy number
                         cncf = box$cncf[,sample_list,drop=FALSE])
  return(box_input_data)
}

runMCMCForABox <- function(box, 
                           n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                           inits = list(".RNG.name" = "base::Wichmann-Hill",
                                        ".RNG.seed" = 123),
                           params = c("z", "w", "ystar"),
                           max_K = 5, model_type = "simple",
                           beta.prior = FALSE,
                           drop_zero = TRUE) {
  # # # returns samps_list 
  # box_input_data <- list(I = nrow(box$y),
  #                        S = ncol(box$y),
  #                        y = box$y,
  #                        n = box$n,
  #                        m = box$m,
  #                        icn = box$tcn, # integer copy number
  #                        cncf = box$cncf)
  
  # modify box_input_data so it only contain non-zero samples
  # if (drop_zero) {
  # sample_list = vector()
  # for (j in 1:ncol(box$y)) {
  #   if (strsplit(box$pattern, "")[[1]][j] == "1") {
  #     sample_list <- append(sample_list, j)
  #   }
  # }
  # # print(temp_box)
  # box_input_data <- list(I = nrow(box$y),
  #                        S = length(sample_list),
  #                        y = box$y[,sample_list,drop=FALSE],
  #                        n = box$n[,sample_list,drop=FALSE],
  #                        m = box$m[,sample_list,drop=FALSE],
  #                        icn = box$tcn[,sample_list,drop=FALSE], # integer copy number
  #                        cncf = box$cncf[,sample_list,drop=FALSE])
  # box_input_data$y <- box$y[,sample_list,drop=FALSE]
  # box_input_data$n <- box$n[,sample_list,drop=FALSE]
  # box_input_data$icn <- box$icn[,sample_list,drop=FALSE]
  # box_input_data$m <- box$m[,sample_list,drop=FALSE]
  # box_input_data$cncf <- box$cncf[,sample_list,drop=FALSE]
  # box_input_data$S <- length(sample_list)
  # }
  box_input_data <- getBoxInputData(box)
  
  extdir <- system.file("extdata", package="pictograph2")
  # if (nrow(box$y) == 1) {
  #   # TO-DO
  #   jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1_I1.jags")
  #   # box_input_data$I <- NULL
  # } else {
  #   # TO-DO
  #   jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1.jags")
  # }
  
  # if (model_type == "simple") {
  #   # TO-DO
  #   jags.file <- file.path(extdir, "model-test.jags") # fixes order of CCFs in one sample, not spike and slab 
  # } else if (model_type == "spike_and_slab") {
  #   # TO-DO
  #   jags.file <- file.path(extdir, "spike_and_slab_purity_ident.jags") # fixing order of CCFs in one sample
  # } else stop("provide model_type either 'spike_and_slab' or 'simple'")
  
  # choose sample in which mutations are present
  sample_to_sort <- which(colSums(box_input_data$y) > 0)[1]
  # box_input_data$sample_to_sort <- sample_to_sort
  
  jags.file.K1 <- file.path(extdir, "spike_and_slab_K1.jags")
  jags.file <- file.path(extdir, "spike_and_slab.jags")
  samps_K1 <- runMCMC(box_input_data, 1, jags.file.K1,
                      inits, params, n.iter=n.iter, thin=thin, n.burn=n.burn)
  
  if(box_input_data$S == 1) {
    colnames(samps_K1[[1]])[which(colnames(samps_K1[[1]]) == "w")] <- "w[1,1]"
  }
  
  # if (drop_zero) {
  samps_K1 <- reverseDrop(samps_K1, box$pattern, n.iter)
  # }
  # 
  # if(box$I == 1) {
  #   colnames(samps_K1[[1]])[which(colnames(samps_K1[[1]]) == "z")] <- "z[1]"
  # }
  
  # Max number of clusters cannot be more than number of mutations
  # max_K <- min(max_K, length(box$mutation_indices)) # already done in the main step
  if (max_K > 1) {
    
    box_input_data$sample_to_sort <- sample_to_sort
    
    samps_2 <- parallel::mclapply(2:max_K,
                                  function(k) runMCMC(box_input_data, k,
                                                      jags.file, inits, params,
                                                      n.iter=n.iter, thin=thin,
                                                      n.burn=n.burn,
                                                      beta.prior=beta.prior),
                                  mc.cores=mc.cores)
    
    # if (drop_zero) {
    for (i in seq_len(length(samps_2))) {
      samps_2[[i]] <- reverseDrop(samps_2[[i]], box$pattern, n.iter)
    }
    # }
    
    samps_list <- c(list(samps_K1), samps_2)
    names(samps_list) <- paste0("K", 1:max_K)
    return(samps_list)
    
  } else {
    names(samps_K1) <- "K1"
    return(samps_K1)
  }
  
}

runMutSetMCMC <- function(temp_box, 
                          n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                          inits = list(".RNG.name" = "base::Wichmann-Hill",
                                       ".RNG.seed" = 123),
                          temp_max_K = 5,
                          model_type = "spike_and_slab",
                          params = c("z", "w", "ystar"),
                          beta.prior = FALSE,
                          drop_zero = FALSE,
                          min_mutation_per_cluster = 1,
                          cluster_diff_thresh=0.05) {
  #beta.prior not used
  # warning("beta prior not used or updated in runMCMC")
  # Run MCMC
  if (temp_max_K == 1) {
    model_type = "simple"
    beta.prior = FALSE
  } 
  
  # WORKING PROGRESS
  temp_samps_list <- runMCMCForABox(temp_box,
                                    n.iter = n.iter, n.burn = n.burn, 
                                    thin = thin, mc.cores = mc.cores,
                                    inits = inits,
                                    params = params,
                                    max_K = temp_max_K,
                                    model = model_type,
                                    beta.prior = beta.prior,
                                    drop_zero = drop_zero)
  
  # Format chains
  if (drop_zero && nrow(temp_box$y) == 1) {
    samps_list <- list(formatChains(temp_samps_list))
    names(samps_list) <- "K1"
  } else {
    samps_list <- parallel::mclapply(temp_samps_list, formatChains,
                                     mc.cores = mc.cores)
  }
  
  # check whether 1) number of mutations per cluster is at least min_mutation_per_cluster 2) difference between any two cluster less than cluster_diff_thresh 
  filtered_samps_list <- filterK(samps_list, min_mutation_per_cluster = min_mutation_per_cluster,
                                 cluster_diff_thresh = cluster_diff_thresh)
  # filtered_samps_list <- samps_list
  # Calculate BIC
  K_tested <- seq_len(length(filtered_samps_list))
  if (temp_max_K > 1) {
    box_indata <- getBoxInputData(temp_box)
    bic_vec <- unname(unlist(parallel::mclapply(filtered_samps_list,
                                                function(chains) calcChainBIC(chains, box_indata, temp_box$pattern),
                                                mc.cores = mc.cores)))
    bic_tb <- tibble(K_tested = K_tested,
                     BIC = bic_vec)
    best_chains <- samps_list[[which.min(bic_vec)]]
    res_list <- list(all_chains = samps_list,
                     BIC = bic_tb,
                     best_chains = best_chains,
                     best_K = which.min(bic_vec))
  } else {
    # only 1 variant, so must be 1 cluster and don't need to check BIC
    res_list <- list(all_chains = filtered_samps_list,
                     BIC = NA,
                     best_chains = filtered_samps_list[[1]],
                     best_K = 1)
  }
  # res_list <- list(all_chains = filtered_samps_list)
  return(res_list)
}

filterK <- function(samps_list, min_mutation_per_cluster=1, cluster_diff_thresh=0.05) {
  filtered_samps_list <- list()
  toBreak = F
  for (i in seq_len(length(samps_list))) {
    k = as.numeric(gsub("\\D", "", names(samps_list)[i]))
    if (k > 1) {
      w_chain = samps_list[[i]]$w_chain
      z_chain = samps_list[[i]]$z_chain
      clusterTable = writeClusterAssignmentsTable(z_chain)
      # check whether all cluster contains at least one mutation
      if (length(unique(clusterTable$Cluster))==k) {
        # check whether all cluster contains at least min_mutation_per_cluster mutations
        if (any(table(clusterTable$Cluster) < min_mutation_per_cluster)) {
          break
        }
      } else {
        break
      }
      mcfTable = writeClusterCCFsTable(w_chain)
      # check whether mcf for any cluster is less than cluster_diff_thresh in all samples
      for (j1 in seq_len(k)) {
        if (all(mcfTable[j1,] < cluster_diff_thresh)) {
          toBreak = T
        }
      }
      # check whether mcf difference between any two clusters less than cluster_diff_thresh in all samples
      for (j1 in seq_len(k-1)) {
        for (j2 in seq(j1+1, k)) {
          diff = abs(mcfTable[j1,] - mcfTable[j2,])
          if (all(diff < cluster_diff_thresh)) {
            toBreak = T
          }
        }
      }
    }
    if (toBreak) { break }
    filtered_samps_list[[names(samps_list)[i]]] <- samps_list[[i]]
  }
  return(filtered_samps_list)
}

formatChains <- function(samps) {
  temp_z <- get.parameter.chain("z", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  temp_w <- get.parameter.chain("w", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  if (nrow(temp_w) == 0) {
    temp_w <- get.parameter.chain("w", ggmcmc::ggs(samps) %>% mutate(Parameter = gsub("w","w[1,1]",Parameter))) %>%
      mutate(Parameter = as.character(Parameter))
  }
  temp_ystar <- get.parameter.chain("ystar", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  samps_list_formatted <- list(w_chain = temp_w,
                               z_chain = temp_z,
                               ystar_chain = temp_ystar)
  return(samps_list_formatted)
}

get.parameter.chain <- function(param, chains) {
  chains[grep(paste0(param, "\\["), chains$Parameter), ]
}
