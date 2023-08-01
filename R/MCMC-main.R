#' run MCMC with added subclonal copy number
#' @export
mcmcMain <- function(max_K = 5, min_mutation_per_cluster=1, iterations=5) {
  # read in files
  data <- importFiles('./inst/extdata/sim_v2_snv.csv', './inst/extdata/sim_v2_cn.csv')
  
  for (iteration in seq_len(iterations)) {
    if (iteration == 1) {
      # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
      cna_init <- estimateCNA(data)
      
      # estimation of copy number cellular fraction
      cncf_init <- estimateCNCF(data, cna_init)
    } else {
      if (all(cna_init==cna_update)) {
        break
      } else {
        cna_init = cna_update
        cncf_init = cncf_update
      }
    }
    
    # assign integer cna and cncf to each mutation
    cna <- data$overlap %*% cna_init + 2*(ifelse(rowSums(data$overlap)==0, 1, 0))
    cncf <- data$overlap %*% cncf_init
    
    # estimate multiplicity
    m_est <- estimateMultiplicity1(data, cna, cncf)
    
    # update data
    input_data <- list(y=data$y,
                       n=data$n,
                       m=m_est,
                       cncf=cncf,
                       tcn=cna)
    
    # 1. separate mutations by sample presence
    sep_list <- separateMutationsBySamplePresence(input_data)
    
    # 2. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
    all_set_results <- runMCMCForAllBoxes(sep_list, max_K = 5, min_mutation_per_cluster = 1, 
                                          n.iter = 5000, n.burn = 1000, thin = 10, mc.cores = 4)
    
    # 3. pick K: most common or min_BIC
    set_k_choices <- writeSetKTable(all_set_results)
    
    # 4. collect best chains
    best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$chosen_K)
    chains <- mergeSetChains(best_set_chains, input_data)
    
    # plotChainsCCF(chains$w_chain)
    # plotCCFViolin(chains$w_chain, chains$z_chain, indata = input_data)
    # plotClusterAssignmentProbVertical(chains$z_chain, chains$w_chain)
    
    # re-estimate cncf by assigning cna to mcf clusters
    cncf_update <- reassignCNCF(cncf_init, chains$w_chain)
    warning("mcmcMain: re-assign cluster mcf after merging CNA")
    cna_update <- reassignCNA(cncf_update, data$tcn)
  }
  
  
  # # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  # cna_init <- estimateCNA(data)
  # 
  # # estimation of copy number cellular fraction
  # cncf_init <- estimateCNCF(data, cna_init)
  # 
  # # assign integer cna and cncf to each mutation
  # cna <- data$overlap %*% cna_init + 2*(ifelse(rowSums(data$overlap)==0, 1, 0))
  # cncf <- data$overlap %*% cncf_init
  # 
  # # estimate multiplicity
  # m_est <- estimateMultiplicity1(data, cna, cncf)
  # 
  # # update data
  # input_data <- list(y=data$y,
  #              n=data$n,
  #              m=m_est,
  #              cncf=cncf,
  #              tcn=cna)
  # 
  # # 1. separate mutations by sample presence
  # sep_list <- separateMutationsBySamplePresence(input_data)
  # 
  # # 2. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
  # all_set_results <- runMCMCForAllBoxes(sep_list, max_K = 5, min_mutation_per_cluster = 1, 
  #                                       n.iter = 5000, n.burn = 1000, thin = 10, mc.cores = 4)
  # 
  # # 3. pick K: most common or min_BIC
  # set_k_choices <- writeSetKTable(all_set_results)
  # 
  # # 4. collect best chains
  # best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$chosen_K)
  # chains <- mergeSetChains(best_set_chains, input_data)
  # 
  # # plotChainsCCF(chains$w_chain)
  # # plotCCFViolin(chains$w_chain, chains$z_chain, indata = input_data)
  # # plotClusterAssignmentProbVertical(chains$z_chain, chains$w_chain)
  # 
  # # re-estimate cncf by assigning cna to mcf clusters
  # cncf_update <- reassignCNCF(cncf_init, chains$w_chain)
  # cna_update <- reassignCNA(cncf_update, data$tcn)
  # 
  writeClusterCCFsTable(chains$w_chain)
  writeClusterAssignmentsTable(chains$z_chain, Mut_ID = data$MutID, cncf=cncf_update)
  
  generateAllTrees(chains$w_chain, lineage_precedence_thresh = 0.1, sum_filter_thresh = 0.2)
  
  scores <- calcTreeScores(chains$w_chain, all_spanning_trees)
  # scores
  # # highest scoring tree
  best_tree <- all_spanning_trees[[which.max(scores)]]
  
  # choose color palette; default is viridis
  color_palette = viridis::viridis
  # # plot tree
  plotTree(best_tree, palette = color_palette)
  plotEnsembleTree(all_spanning_trees, palette = color_palette)
  
  subclone_props <- calcSubcloneProportions(w_mat, best_tree)
  plotSubclonePie(subclone_props, palette=color_palette, sample_names=colnames(input_data$y))
  plotSubcloneBar(subclone_props, palette=color_palette, sample_names=colnames(input_data$y))
  
}

