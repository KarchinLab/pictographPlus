#' run MCMC with added subclonal copy number
#' 
#' @export
mcmcMain <- function(mutation_file,
                     copy_number_file,
                     outputDir,
                     SNV_file=NULL,
                     stat_file=NULL, 
                     cytoband_file=NULL, 
                     pval=0.05,
                     sim_iter=100,
                     max_K = 3, 
                     min_mutation_per_cluster=5, 
                     iterations=5, 
                     min_mutation_per_box=10, 
                     n.iter=5000, 
                     n.burn=1000, 
                     thin=10, 
                     mc.cores=8, 
                     model_type="spike_and_slab", 
                     beta.prior=FALSE, 
                     drop_zero=TRUE, 
                     inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
                     cluster_diff_thresh=0.05,
                     alt_reads_thresh = 0, 
                     vaf_thresh = 0, 
                     cnv_max_dist=2000, 
                     cnv_max_percent=0.30, 
                     tcn_normal_range=c(1.8, 2.2), 
                     smooth_cnv=T, 
                     autosome=T) {
  
  data <- importFiles(mutation_file, 
                      copy_number_file, 
                      outputDir, 
                      SNV_file=SNV_file, 
                      stat_file=stat_file, 
                      cytoband_file=cytoband_file, 
                      alt_reads_thresh=alt_reads_thresh, 
                      vaf_thresh=vaf_thresh, 
                      cnv_max_dist=cnv_max_dist, 
                      cnv_max_percent=cnv_max_percent, 
                      tcn_normal_range=tcn_normal_range, 
                      smooth_cnv=smooth_cnv, 
                      autosome=autosome, 
                      mc.cores=mc.cores, 
                      pval=pval,
                      sim_iter=sim_iter)


  # # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  # cna_init <- estimateCNA(data)
  #   
  # # estimation of copy number cellular fraction
  # cncf_init <- estimateCNCF(data, cna_init)
  # 
  # 
  # # assign integer cna and cncf to each mutation
  # icn <- data$overlap %*% cna_init + 2*(ifelse(rowSums(data$overlap)==0, 1, 0))
  # cncf <- data$overlap %*% cncf_init
  # 
  # # estimate multiplicity
  # m_est <- estimateMultiplicity1(data, icn, cncf)
  
  # get index of SNVs in CN-neutral region (no LOH) or CNA events
  index = which(rowSums(data$tcn)==0 | data$is_cn==1)
  # update data
  input_data <- list(y=data$y[index,],
                     n=data$n[index,],
                     tcn=data$tcn[index,],
                     is_cn=data$is_cn[index],
                     MutID=data$MutID[index])
  
  # 1. separate mutations by sample presence
  sep_list <- separateMutationsBySamplePresence(input_data, min_mutation_per_box)
  
  # 2. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
  all_set_results <- runMCMCForAllBoxes(sep_list, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                        n.iter = 5000, n.burn = 1000, thin = 10, mc.cores = 8, model_type = "type1")
  
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
  
  
  
  
  
  
  generateAllTrees(chains$w_chain, lineage_precedence_thresh = 0.02, sum_filter_thresh = 0.1)
  ccfTable = writeClusterCCFsTable(chains$w_chain)
  
  # write.table(ccfTable, file=paste(outputDir, "ccf.csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  
  clusterAssingmentTable = writeClusterAssignmentsTable(chains$z_chain, w_chain=chains$w_chain, Mut_ID = data$MutID, cncf=cncf_update)
  
  for (i in seq_len(nrow(clusterAssingmentTable))) {
    # rename chromosome name if CNA in a cluster
    if (clusterAssingmentTable[i,]$Mut_ID %in% rownames(data$tcn)) {
      change = "neutral"
      if (mean(cna_update[(clusterAssingmentTable[i,]$Mut_ID),]) < 2) {
        change = "del"
      }
      if (mean(cna_update[(clusterAssingmentTable[i,]$Mut_ID),]) > 2) {
        change = "dup"
      }
      
      ######################################
      # Add code to add cytoband information
      ######################################
      new_name = paste(change, clusterAssingmentTable[i,]$Mut_ID, sep = ":")
      clusterAssingmentTable[i,]$Mut_ID = new_name
    }
  }
  
  # write.table(clusterAssingmentTable, file=paste(outputDir, "clusterAssign.csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  
  scores <- calcTreeScores(chains$w_chain, all_spanning_trees)

  # highest scoring tree
  best_tree <- all_spanning_trees[[which.max(scores)]]
  best_tree <- all_spanning_trees[[12]] # for HTAN-IPMN MCL111_001

  # plot tree
  plotTree(best_tree, palette = viridis::viridis)
  # plotEnsembleTree(all_spanning_trees, palette = color_palette)
  
  subclone_props <- calcSubcloneProportions(w_mat, best_tree)
  plotSubclonePie(subclone_props, sample_names=colnames(input_data$y))
  plotSubcloneBar(subclone_props, sample_names=colnames(input_data$y))
  # save.image(file=paste(outputDir, "PICTograph2.RData", sep=""))
}

