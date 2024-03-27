#' run MCMC with added subclonal copy number
#' 
#' @export
mcmcMain <- function(mutation_file,
                     copy_number_file,
                     outputDir,
                     SNV_file=NULL,
                     stat_file=NULL, 
                     cytoband_file=NULL, 
                     sample_presence=FALSE,
                     ploidy=2,
                     pval=0.05,
                     max_K = 5, 
                     min_mutation_per_cluster=5, 
                     n.iter=5000, 
                     n.burn=1000, 
                     thin=10, 
                     mc.cores=8, 
                     inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
                     cluster_diff_thresh=0.05,
                     alt_reads_thresh = 0, 
                     vaf_thresh = 0, 
                     cnv_max_dist=2000, 
                     cnv_max_percent=0.30, 
                     tcn_normal_range=c(1.8, 2.2), 
                     smooth_cnv=F, 
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
                      pval=pval)
  
  data <- assign("data", data, envir = .GlobalEnv)
  
  min_mutation_per_cluster <- max(floor(nrow(data$y)/50), min_mutation_per_cluster)
  
  if (ncol(data$y)==1) {
    min_mutation_per_cluster <- max(floor(nrow(data$y)/25), min_mutation_per_cluster)
  } 

  if (sample_presence) {
    
    message("Using sample presence")
    ##############################################################################
    #              MCMC 1: SSMs (CN-neutral region) and all CNAs                 #
    ##############################################################################
    # 1. get index of SNVs in CN-neutral region (no LOH) or CNA events
    index = which(rowSums(data$tcn)==0 | data$is_cn==1)
    # 2. update input_data
    input_data <- list(y=data$y[index,,drop=FALSE],
                       n=data$n[index,,drop=FALSE],
                       tcn=data$tcn[index,,drop=FALSE],
                       is_cn=data$is_cn[index],
                       MutID=data$MutID[index])
    
    # 3. separate mutations by sample presence
    sep_list <- separateMutationsBySamplePresence(input_data)
    
    # 4. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
    #    Note: model_type = "type1"
    all_set_results <- runMCMCForAllBoxes(sep_list, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                          cluster_diff_thresh = cluster_diff_thresh, inits = inits,
                                          n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type1")
    

    # 5. pick K: most common or min_BIC
    set_k_choices <- writeSetKTable(all_set_results)
    
    # 6. collect best chains
    best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$chosen_K)
    chains <- mergeSetChains(best_set_chains, input_data)
    
    # 7. find integer copy number and multiplicity
    data$mtp <- findM(data, input_data, chains)
    data$icn <- findIcn(data, input_data, chains)
    data$cncf <- findCncf(data, input_data, chains)
    
    ##############################################################################
    #              MCMC 2: all SSMs and all CNAs                 #
    ##############################################################################
    
    # 8. collect data for the second chain
    input_data <- list(y=data$y,
                       n=data$n,
                       tcn=data$tcn,
                       is_cn=data$is_cn,
                       mtp=data$mtp,
                       icn=data$icn,
                       cncf=data$cncf,
                       MutID=data$MutID)
    
    input_data <- assign("input_data", input_data, envir = .GlobalEnv)
    # 9. separate mutations by sample presence
    sep_list <- separateMutationsBySamplePresence(input_data)
    
    # 10. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
    #    Note: model_type = "type1"
    all_set_results <- runMCMCForAllBoxes(sep_list, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                          cluster_diff_thresh = cluster_diff_thresh, inits = inits,
                                          n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
    
    all_set_results <- assign("all_set_results", all_set_results, envir = .GlobalEnv)
    
    # 11. pick K: most common or min_BIC
    set_k_choices <- writeSetKTable(all_set_results)
    set_k_choices <- assign("set_k_choices", set_k_choices, envir = .GlobalEnv)
    
    # 12. collect best chains
    best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$chosen_K)
    chains <- mergeSetChains(best_set_chains, input_data)
  
  } else {
    ##############################################################################
    #             MCMC chain for all; no sample presence                         #
    ##############################################################################
    max_K <- max_K * ncol(data$y)
    
    input_data <- list(y=data$y,
                       n=data$n,
                       tcn=data$tcn,
                       is_cn=data$is_cn,
                       q=data$q,
                       MutID=data$MutID)
    
    all_set_results <- runMCMCForAllBoxes(input_data, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                          cluster_diff_thresh = cluster_diff_thresh, inits = inits,
                                          n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type3")
    
    all_set_results <- assign("all_set_results", all_set_results, envir = .GlobalEnv)
    
    # 1. pick K: chosen K is min_BIC
    set_k_choices <- writeSetKTable(all_set_results)
    set_k_choices <- assign("set_k_choices", set_k_choices, envir = .GlobalEnv)
    
    # 2. collect best chains
    best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$chosen_K)

    chains <- mergeSetChains(best_set_chains, input_data)
    
    
  }
  
  
  png(paste(outputDir, "mcf.png", sep="/"))
  plotChainsMCF(chains$mcf_chain)
  dev.off()
  # plotMCFViolin(chains$mcf_chain, chains$z_chain, indata = input_data)
  # plotClusterAssignmentProbVertical(chains$z_chain, chains$mcf_chain)
  
  mcfTable = writeClusterMCFsTable(chains$mcf_chain)
  colnames(mcfTable)=c("Cluster",c(colnames(data$y)))
  # mcfTable
  
  write.table(mcfTable, file=paste(outputDir, "mcf.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  clusterAssingmentTable = writeClusterAssignmentsTable(chains$z_chain, Mut_ID = input_data$MutID)
  # clusterAssingmentTable
  
  icnTable <- writeIcnTable(chains$icn_chain, Mut_ID = input_data$MutID)
  write.table(icnTable, file=paste(outputDir, "icn_all.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
   
  multiplicityTable <- writeMultiplicityTable(chains$m_chain, Mut_ID = input_data$MutID)
  write.table(multiplicityTable, file=paste(outputDir, "multiplicity_all.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  icnTableCN <- icnTable[data$is_cn==1,]
  multiplicityTableCN <- multiplicityTable[data$is_cn==1,]
  icnTableCN$Multiplicity <- multiplicityTableCN$Multiplicity
  write.table(icnTableCN, file=paste(outputDir, "CN_results.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  ### Clean copy number segments by removing segments with icn of 2 and multiplicity of 1
  
  toKeepIndex = c()
  for (i in seq_len(nrow(clusterAssingmentTable))) {
    # rename chromosome name if CNA in a cluster
    
    if (clusterAssingmentTable[i,]$Mut_ID %in% icnTableCN$Mut_ID) {
      # print(clusterAssingmentTable[i,]$Mut_ID)
      
      icnInfo <- icnTableCN %>% filter(Mut_ID==clusterAssingmentTable[i,]$Mut_ID)
      if (icnInfo$icn==2 & icnInfo$Multiplicity==1) {
        # print(icnInfo)
        # stop("stop")
        toKeepIndex <- c(toKeepIndex, 0)
      } else {
        major_cn = max(icnInfo$Multiplicity, icnInfo$icn-icnInfo$Multiplicity)
        clusterAssingmentTable[i,]$Mut_ID <- paste(clusterAssingmentTable[i,]$Mut_ID, ";icn:", icnInfo$icn, ";", "major_cn:", major_cn, sep="")
        # print(change)
        toKeepIndex <- c(toKeepIndex, 1)
        # stop()
      }
      
    } else {
      toKeepIndex <- c(toKeepIndex, 1)
    }
  }
  clusterAssingmentTable$idx <- toKeepIndex
  clusterAssingmentTable <- clusterAssingmentTable %>% filter(toKeepIndex==1) %>% select(-idx)
  write.table(clusterAssingmentTable, file=paste(outputDir, "clusterAssign.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  threshes <- allThreshes()
  # print(threshes)
  
  for (thresh in threshes) {
    generateAllTrees(chains$mcf_chain, lineage_precedence_thresh = thresh[1], sum_filter_thresh = thresh[2])
    if (length(all_spanning_trees) > 0) {
      break
    }
  }

  scores <- calcTreeScores(chains$mcf_chain, all_spanning_trees)

  # # plot all tree with best scores
  # for (i in seq_len(length(which(scores == max(scores))))) {
  #   idx = which(scores == max(scores))[i]
  #   print(idx)
  #   best_tree <- all_spanning_trees[[idx]]
  #   write.table(best_tree, file=paste(outputDir, "/tree", i, ".csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  #   
  #   png(paste(outputDir, "/tree", i, ".png", sep=""))
  #   # plot tree
  #   plotTree(best_tree, palette = viridis::viridis)
  #   # plotEnsembleTree(all_spanning_trees, palette = viridis::viridis)
  #   dev.off()
  #   
  #   cc <- best_tree %>% filter(parent=="root") %>% select(child)
  #   purity <- mcfTable %>% filter(Cluster %in% cc$child) %>% summarise(across(everything(), sum)) %>% select(-Cluster)
  #   colnames(purity) <- colnames(data$y)
  #   write.table(purity, file=paste(outputDir, "/purity", i, ".csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  #   
  #   subclone_props <- calcSubcloneProportions(mcf_mat, best_tree)
  #   rownames(subclone_props) = mcfTable$Cluster
  #   colnames(subclone_props) = colnames(data$y)
  #   
  #   write.csv(subclone_props, file=paste(outputDir, "/subclone_proportion", i, ".csv", sep=""), quote = FALSE)
  #   
  #   png(paste(outputDir, "/subclone_props", i, ".png", sep=""))
  #   plotSubclonePie(subclone_props, sample_names=colnames(input_data$y))
  #   dev.off()
  # }

  
  # highest scoring tree
  best_tree <- all_spanning_trees[[which(scores == max(scores))[length(which(scores == max(scores)))]]]
  write.table(best_tree, file=paste(outputDir, "tree.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)

  # plot tree
  png(paste(outputDir, "tree.png", sep="/"))

  plotTree(best_tree, palette = viridis::viridis)
  # plotEnsembleTree(all_spanning_trees, palette = viridis::viridis)
  dev.off()

  cc <- best_tree %>% filter(parent=="root") %>% select(child)
  purity <- mcfTable %>% filter(Cluster %in% cc$child) %>% summarise(across(everything(), sum)) %>% select(-Cluster)
  colnames(purity) <- colnames(data$y)
  write.table(purity, file=paste(outputDir, "purity.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)

  subclone_props <- calcSubcloneProportions(mcf_mat, best_tree)
  rownames(subclone_props) = mcfTable$Cluster
  colnames(subclone_props) = colnames(data$y)

  write.csv(subclone_props, file=paste(outputDir, "subclone_proportion.csv", sep="/"), quote = FALSE)

  png(paste(outputDir, "subclone_props.png", sep="/"))
  plotSubclonePie(subclone_props, sample_names=colnames(input_data$y))
  dev.off()
  # plotSubclonePie(subclone_props, sample_names=colnames(input_data$y))
  # plotSubcloneBar(subclone_props, sample_names=colnames(input_data$y))
  
  save.image(file=paste(outputDir, "PICTograph2.RData", sep="/"))
}

findM <- function(data, input_data, chains) {
  mTable <- writeMultiplicityTable(chains$m_chain, Mut_ID = input_data$MutID)
  mTable <- mTable %>% 
    mutate(order_idx = match(mTable$Mut_ID, input_data$MutID)) %>% 
    arrange(order_idx) %>% 
    select(-order_idx)
  mTable <- mTable[which(input_data$is_cn==1),] 
  mTable <- mTable %>%
    mutate(order_idx = match(mTable$Mut_ID, data$MutID[which(data$is_cn==1)])) %>% 
    arrange(order_idx) %>% 
    select(-order_idx)
  
  ifelse(rowSums(data$tcn)==0, 1, data$overlap %*% as.matrix(mTable$Multiplicity))
}

findIcn <- function(data, input_data, chains) {
  icnTable <- writeIcnTable(chains$icn_chain, Mut_ID = input_data$MutID)
  icnTable <- icnTable %>% 
    mutate(order_idx = match(icnTable$Mut_ID, input_data$MutID)) %>% 
    arrange(order_idx) %>% 
    select(-order_idx)
  icnTable <- icnTable[which(input_data$is_cn==1),]
  icnTable <- icnTable %>%
    mutate(order_idx = match(icnTable$Mut_ID, data$MutID[which(data$is_cn==1)])) %>% 
    arrange(order_idx) %>% 
    select(-order_idx)
  
  ifelse(rowSums(data$tcn)==0, 2, data$overlap %*% as.matrix(icnTable$icn))
}

findCncf <- function(data, input_data, chains) {
  cTable = writeClusterAssignmentsTable(chains$z_chain, Mut_ID = input_data$MutID)
  cTable <- cTable %>% 
    mutate(order_idx = match(cTable$Mut_ID, input_data$MutID)) %>% 
    arrange(order_idx) %>% 
    select(-order_idx)
  mcfTable = writeClusterMCFsTable(chains$mcf_chain)
  colnames(mcfTable) <- c("Cluster", colnames(input_data$y))
  tmp <- cTable %>% inner_join(mcfTable, by="Cluster") %>% select(-Cluster)
  tmp1 <- tmp[which(input_data$is_cn==1),]
  tmp1 <- tmp1 %>% 
    mutate(order_idx = match(tmp1$Mut_ID, colnames(data$overlap))) %>% 
    arrange(order_idx) %>% 
    select(-order_idx)
  
  overlap <- data$overlap
  result_tibble <- tibble(row=1:nrow(overlap))
  result_tibble$Mut_ID <- apply(overlap, 1, function(row) {
    cols <- which(row == 1)
    paste(colnames(overlap)[cols], collapse=", ")
  })
  result_tibble$row <- rownames(overlap)
  tmp2 <-tmp1 %>% 
    full_join(result_tibble, by="Mut_ID") %>%
    select(-Mut_ID) %>%
    rename(Mut_ID = row) %>% 
    replace(is.na(.), 0) 
  tmp2 <- tmp2 %>%
    mutate(order_idx = match(tmp2$Mut_ID, data$MutID)) %>% 
    arrange(order_idx) %>% 
    select(-order_idx) 
  
  row_name <- tmp2$Mut_ID
  tmp2 <- tmp2 %>% select(-Mut_ID)
  tmp2 <- as.matrix(tmp2)
  rownames(tmp2) <- row_name
  tmp2
}

allThreshes <- function() {
  threshes <- list() 
  threshes[[1]] <- c(0,0)
  threshes[[2]] <- c(0.1,0)
  threshes[[3]] <- c(0,0.1)
  threshes[[4]] <- c(0.1,0.1)
  threshes[[5]] <- c(0.1,0.2)
  threshes[[6]] <- c(0.2,0.1)
  threshes[[7]] <- c(0.2,0.2)
  threshes[[8]] <- c(0.1,0.3)
  threshes[[9]] <- c(0.3,0.1)
  threshes[[10]] <- c(0.2,0.3)
  threshes[[11]] <- c(0.3,0.2)
  threshes[[12]] <- c(0.3,0.3)
  threshes[[13]] <- c(0.4,0.4)
  threshes
}

# findOverlap <- function(data) {
#   overlap <- data$overlap
#   colnames(overlap) <- sapply(colnames(data$overlap), function(col) {which(rownames(data$overlap)==col)})
#   tmp <- vector("numeric", nrow(overlap))
#   for (i in 1:nrow(overlap)) {
#     tmp[i] <- ifelse(length(which(overlap[i,] == 1)) > 0, as.numeric(names(which(overlap[i,] == 1))[1]),i)
#   }
#   tmp
# }
