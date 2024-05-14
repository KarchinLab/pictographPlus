#' run PICTograph2 in an automated pipeline
#' 
#' run MCMC chains to  infer the clonal evolution of tumors from single or multi-region sequencing data. 
#' This function automatically runs a pipeline of the tool. It models uncertainty of mutation cellular 
#' fraction (MCF) in small somatic mutations (SSMs) and copy number alterations (CNAs), assigning SSMs 
#' and CNAs to subclones using a Bayesian hierarchical model, and reconstruct tumor evolutionary trees 
#' that are constrained based on principles of lineage precedence, sum condition, and optionally by 
#' sample-presence. 
#' 
#' @param mutation_file a csv file that include information for SSMs.
#' @param copy_number_file a csv file that include information for CNA.
#' @param SNV_file a csv file that include information for germline heterozygous SNVs.
#' @param outputDir output directory for saving all files.
#' @param sample_presence whether to use sample presence to separate the mutations. Not applicable if dual_model is set to FALSE and a copy number file is provided.
#' @param score scoring function to estimate the number of clusters. silhouette or BIC.
#' @param max_K user defined maximum number of clusters.
#' @param min_mutation_per_cluster minumum number of mutations in each cluster.
#' @param n.iter number of iterations by JAGS.
#' @param n.burn number of burns by JAGS.
#' @param thin number of thin by JAGS.
#' @param mc.cores number of cores to use for parallel computing; not applicable to windows.
#' @param inits additional parameters by JAGS.
#' @param cluster_diff_thresh threshold to merge two clusters.
#' @export
mcmcMain <- function(mutation_file,
                     copy_number_file=NULL,
                     SNV_file=NULL,
                     outputDir=NULL,
                     sample_presence=TRUE,
                     dual_model=TRUE, # placeholder; dual_model=FALSE still require testing
                     score="silhouette", # either BIC or silhouette
                     ploidy=2, # placeholder
                     pval=0.05, # placeholder
                     max_K = 10, 
                     min_mutation_per_cluster=5, 
                     cluster_diff_thresh=0.05,
                     n.iter=5000, 
                     n.burn=1000, 
                     thin=10, 
                     mc.cores=8, 
                     inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
                     alt_reads_thresh = 0, # placeholder
                     vaf_thresh = 0, # placeholder
                     cnv_max_dist=2000, # placeholder
                     cnv_max_percent=0.30, # placeholder
                     tcn_normal_range=c(1.7, 2.3), # placeholder
                     smooth_cnv=F, # placeholder
                     autosome=T # placeholder
                     ) {
  
  data <- importFiles(mutation_file, 
                      copy_number_file, 
                      outputDir, 
                      SNV_file=SNV_file, 
                      alt_reads_thresh=alt_reads_thresh, 
                      vaf_thresh=vaf_thresh, 
                      cnv_max_dist=cnv_max_dist, 
                      cnv_max_percent=cnv_max_percent, 
                      tcn_normal_range=tcn_normal_range, 
                      smooth_cnv=smooth_cnv, 
                      autosome=autosome, 
                      pval=pval)
  
  # use working directory to save outputs if outputDir is not provided
  if (is.null(outputDir)) {
    outputDir = getwd()
  }
  
  # save upset plot if more than one sample
  if (ncol(data$y) > 1) {
    data_matrix <- ifelse(data$y[data$is_cn==0,]>0, 1, 0)
    png(paste(outputDir, "upsetR.png", sep="/"), res=100)
    print(upset(as.data.frame(data_matrix), text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5), keep.order = T, sets = rev(colnames(data_matrix))))
    dev.off()
  }
  
  data <- assign("data", data, envir = .GlobalEnv)
  
  # use model 2 if only mutation file is provided
  if (is.null(copy_number_file)) { 
    
    if (sample_presence) {
      message("Using sample presence; SSM only")
      input_data <- list(y=data$y,
                         n=data$n,
                         tcn=data$tcn,
                         is_cn=data$is_cn,
                         mtp=data$mtp,
                         icn=data$icn,
                         cncf=data$cncf,
                         MutID=data$MutID,
                         purity=data$purity)
      
      input_data <- assign("input_data", input_data, envir = .GlobalEnv)
      
      # separate mutations by sample presence
      sep_list <- separateMutationsBySamplePresence(input_data)
      
      # For each presence set, run clustering MCMC, calculate silhouette and BIC and choose best K
      all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                            cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                            n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
    } else {
      message("Not using sample presence; SSM only")
      input_data <- list(y=data$y,
                   n=data$n,
                   tcn=data$tcn,
                   is_cn=data$is_cn,
                   mtp=data$mtp,
                   icn=data$icn,
                   cncf=data$cncf,
                   MutID=data$MutID,
                   purity=data$purity)

      input_data <- assign("input_data", input_data, envir = .GlobalEnv)
      
      all_set_results <- runMCMCForAllBoxes(input_data, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                            cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                            n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
    }
  } else {
    # using two step modeling
    if (dual_model) {
      if (sample_presence) {
        
        message("Using sample presence; SSM and CNA")
        ##############################################################################
        #              MCMC 1: SSMs (CN-neutral region) and all CNAs                 #
        ##############################################################################
        # get index of SNVs in CN-neutral region (no LOH) or CNA events
        index = which(rowSums(data$tcn)==0 | data$is_cn==1)
        
        # update input_data
        input_data <- list(y=data$y[index,,drop=FALSE],
                           n=data$n[index,,drop=FALSE],
                           tcn=data$tcn[index,,drop=FALSE],
                           is_cn=data$is_cn[index],
                           MutID=data$MutID[index],
                           purity=data$purity)
        
        # separate mutations by sample presence
        sep_list <- separateMutationsBySamplePresence(input_data)
        
        # For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
        all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type1")
        
        
        # pick K
        set_k_choices <- writeSetKTable(all_set_results)
        
        # collect best chains
        if (score == "silhouette") {
          best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$silhouette_K)
        } else {
          best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$BIC_K)
        }
        
        chains <- mergeSetChains(best_set_chains, input_data)
        
        # 7. find integer copy number and multiplicity
        data$mtp <- findM(data, input_data, chains)
        data$icn <- findIcn(data, input_data, chains)
        data$cncf <- findCncf(data, input_data, chains)
        
        ##############################################################################
        #              MCMC 2: all SSMs and all CNAs                 #
        ##############################################################################
        # collect data for the second model
        input_data <- list(y=data$y,
                           n=data$n,
                           tcn=data$tcn,
                           is_cn=data$is_cn,
                           mtp=data$mtp,
                           icn=data$icn,
                           cncf=data$cncf,
                           MutID=data$MutID,
                           purity=data$purity)
        
        input_data <- assign("input_data", input_data, envir = .GlobalEnv)
        
        # separate mutations by sample presence
        sep_list <- separateMutationsBySamplePresence(input_data)
        
        # For each presence set, run clustering MCMC
        all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
        
      } else {
        
        warning("Not using sample presence; SSM and CNA; NOT TESTED YET")
        ##############################################################################
        #              MCMC 1: SSMs (CN-neutral region) and all CNAs                 #
        ##############################################################################
        # get index of SNVs in CN-neutral region (no LOH) or CNA events
        index = which(rowSums(data$tcn)==0 | data$is_cn==1)
        # update input_data
        input_data <- list(y=data$y[index,,drop=FALSE],
                           n=data$n[index,,drop=FALSE],
                           tcn=data$tcn[index,,drop=FALSE],
                           is_cn=data$is_cn[index],
                           MutID=data$MutID[index],
                           purity=data$purity)
        
        # For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
        all_set_results <- runMCMCForAllBoxes(input_data, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type1")
        # pick K
        set_k_choices <- writeSetKTable(all_set_results)
        
        # collect best chains
        if (score == "silhouette") {
          best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$silhouette_K)
        } else {
          best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$BIC_K)
        }
        chains <- mergeSetChains(best_set_chains, input_data)
        
        # find integer copy number and multiplicity
        data$mtp <- findM(data, input_data, chains)
        data$icn <- findIcn(data, input_data, chains)
        data$cncf <- findCncf(data, input_data, chains)
        
        ##############################################################################
        #              MCMC 2: all SSMs and all CNAs                 #
        ##############################################################################
        
        # collect data for the second chain
        input_data <- list(y=data$y,
                           n=data$n,
                           tcn=data$tcn,
                           is_cn=data$is_cn,
                           mtp=data$mtp,
                           icn=data$icn,
                           cncf=data$cncf,
                           MutID=data$MutID,
                           purity=data$purity)
        
        input_data <- assign("input_data", input_data, envir = .GlobalEnv)

        all_set_results <- runMCMCForAllBoxes(input_data, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
        
      }
      
    } else {
      message("Not using sample presence; Using single model for both SSM and CNA")
      ##############################################################################
      #             MCMC chain for all; no sample presence                         #
      ##############################################################################
      input_data <- list(y=data$y,
                         n=data$n,
                         tcn=data$tcn,
                         is_cn=data$is_cn,
                         q=data$q,
                         MutID=data$MutID,
                         purity=data$purity)
      
      input_data <- assign("input_data", input_data, envir = .GlobalEnv)
      
      all_set_results <- runMCMCForAllBoxes(input_data, sample_presence=FALSE, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                            cluster_diff_thresh = cluster_diff_thresh, inits = inits,
                                            n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type3")
    }
  }
  
  all_set_results <- assign("all_set_results", all_set_results, envir = .GlobalEnv)
  
  # pick K: silhouette or BIC
  set_k_choices <- writeSetKTable(all_set_results)
  set_k_choices <- assign("set_k_choices", set_k_choices, envir = .GlobalEnv)
  
  # collect best chains
  if (score=="silhouette") {
    best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$silhouette_K)
  } else {
    best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$BIC_K)
  }
  chains <- mergeSetChains(best_set_chains, input_data)
  
  # plot MCMC tracing 
  png(paste(outputDir, "mcf.png", sep="/"))
  print(
    plotChainsMCF(chains$mcf_chain)
  )
  dev.off()
  
  # plot violin plot
  # png(paste(outputDir, "violin.png", sep="/"))
  # print(
  #   plotMCFViolin(chains$mcf_chain, chains$z_chain, indata = input_data)
  # )
  # dev.off()
  
  # write mcf table
  mcfTable = writeClusterMCFsTable(chains$mcf_chain)
  colnames(mcfTable)=c("Cluster",c(colnames(data$y)))
  write.table(mcfTable, file=paste(outputDir, "mcf.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # write cluster assignment table
  clusterAssingmentTable = writeClusterAssignmentsTable(chains$z_chain, Mut_ID = input_data$MutID)
  
  # record estimated icn and multiplicity information
  icnTable <- writeIcnTable(chains$icn_chain, Mut_ID = input_data$MutID)
  write.table(icnTable, file=paste(outputDir, "icn_all.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  multiplicityTable <- writeMultiplicityTable(chains$m_chain, Mut_ID = input_data$MutID)
  write.table(multiplicityTable, file=paste(outputDir, "multiplicity_all.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  icnTableCN <- icnTable[data$is_cn==1,]
  multiplicityTableCN <- multiplicityTable[data$is_cn==1,]
  icnTableCN$Multiplicity <- multiplicityTableCN$Multiplicity
  write.table(icnTableCN, file=paste(outputDir, "CN_results.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # generate trees using different set of thresholds until at least one tree is available
  threshes <- allThreshes()
  for (thresh in threshes) {
    generateAllTrees(chains$mcf_chain, data$purity, lineage_precedence_thresh = thresh[1], sum_filter_thresh = thresh[2])
    if (length(all_spanning_trees) > 0) {
      break
    }
  }

  if (is.null(copy_number_file)) {
    cncfTable <- data$cncf
  } else {
    cncfTable <- findCncf(data, input_data, chains)
  }
  # scores <- calcTreeScores(chains$mcf_chain, all_spanning_trees, purity=data$purity)
  scores <- calculateTreeScoreMutations(chains$mcf_chain, data, icnTable, cncfTable, multiplicityTable, clusterAssingmentTable, data$purity, all_spanning_trees)
  scores <- assign("scores", scores, envir = .GlobalEnv)
  
  # plot all possible trees
  plotAllTrees(outputDir, scores, all_spanning_trees, mcfTable, data)
  
  # highest scoring tree
  best_tree <- all_spanning_trees[[which(scores == max(scores))[length(which(scores == max(scores)))]]]
  write.table(best_tree, file=paste(outputDir, "tree.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # plot best and ensemble tree
  if (nrow(best_tree) >1) {
    png(paste(outputDir, "tree.png", sep="/"))
    plotTree(best_tree, palette = viridis::viridis)
    dev.off()
  }
  
  if (length(which(scores == max(scores))) > 1) {
    png(paste(outputDir, "tree_ensemble.png", sep="/"))
    plotEnsembleTree(all_spanning_trees[which(scores == max(scores))], palette = viridis::viridis)
    dev.off()
  }
  
  # estimate purity
  cc <- best_tree %>% filter(parent=="root") %>% select(child)
  purity <- mcfTable %>% filter(Cluster %in% cc$child) %>% summarise(across(everything(), sum)) %>% select(-Cluster)
  colnames(purity) <- colnames(data$y)
  write.table(purity, file=paste(outputDir, "purity.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # estimate subclone proportion
  subclone_props <- calcSubcloneProportions(mcf_mat, best_tree)
  rownames(subclone_props) = mcfTable$Cluster
  colnames(subclone_props) = colnames(data$y)

  write.csv(subclone_props, file=paste(outputDir, "subclone_proportion.csv", sep="/"), quote = FALSE)

  png(paste(outputDir, "subclone_props.png", sep="/"))
  print(plotSubclonePie(subclone_props, sample_names=colnames(input_data$y)))
  dev.off()

  # Clean copy number segments by removing segments with icn of 2 and multiplicity of 1
  toKeepIndex = c()
  for (i in seq_len(nrow(clusterAssingmentTable))) {
    if (clusterAssingmentTable[i,]$Mut_ID %in% icnTableCN$Mut_ID) {
      icnInfo <- icnTableCN %>% filter(Mut_ID==clusterAssingmentTable[i,]$Mut_ID)
      if (icnInfo$icn==2 & icnInfo$Multiplicity==1) {
        toKeepIndex <- c(toKeepIndex, 0)
      } else {
        major_cn = max(icnInfo$Multiplicity, icnInfo$icn-icnInfo$Multiplicity)
        clusterAssingmentTable[i,]$Mut_ID <- paste(clusterAssingmentTable[i,]$Mut_ID, ";icn:", icnInfo$icn, ";", "major_cn:", major_cn, sep="")
        toKeepIndex <- c(toKeepIndex, 1)
      }
    } else {
      toKeepIndex <- c(toKeepIndex, 1)
    }
  }
  clusterAssingmentTable$idx <- toKeepIndex
  clusterAssingmentTable <- clusterAssingmentTable %>% filter(toKeepIndex==1) %>% select(-idx)
  write.table(clusterAssingmentTable, file=paste(outputDir, "clusterAssign.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # save all data
  save.image(file=paste(outputDir, "PICTograph2.RData", sep="/"))
  
}


#' Plot all trees with the highest scores
plotAllTrees <- function(outputDir, scores, all_spanning_trees, mcfTable, data) {
  # plot all tree with best scores
  
  outputDir = paste(outputDir, "all_trees", sep = "/")
  suppressWarnings(dir.create(outputDir))
  
  for (i in seq_len(length(which(scores == max(scores))))) {
    idx = which(scores == max(scores))[i]
    
    best_tree <- all_spanning_trees[[idx]]
    if (nrow(best_tree) >1 ) {
      write.table(best_tree, file=paste(outputDir, "/tree", i, ".csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  
      png(paste(outputDir, "/tree", i, ".png", sep=""))
      # plot tree
      plotTree(best_tree, palette = viridis::viridis)
      # plotEnsembleTree(all_spanning_trees, palette = viridis::viridis)
      dev.off()
  
      cc <- best_tree %>% filter(parent=="root") %>% select(child)
      purity <- mcfTable %>% filter(Cluster %in% cc$child) %>% summarise(across(everything(), sum)) %>% select(-Cluster)
      colnames(purity) <- colnames(data$y)
      write.table(purity, file=paste(outputDir, "/tree_", i, "_purity.csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  
      subclone_props <- calcSubcloneProportions(mcf_mat, best_tree)
      rownames(subclone_props) = mcfTable$Cluster
      colnames(subclone_props) = colnames(data$y)
  
      write.csv(subclone_props, file=paste(outputDir, "/tree_", i, "_subclone_proportion.csv", sep=""), quote = FALSE)
  
      png(paste(outputDir, "/tree_", i, "_subclone_proportion.png", sep=""))
      print(plotSubclonePie(subclone_props, sample_names=colnames(input_data$y)))
      dev.off()
    }
  }
}

#' find multiplicity for each mutation
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

#' find integer copy number for each mutaiton
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

#' find cncf for each mutation
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

#'  defines the thresholds to be used for tree building
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
  threshes[[13]] <- c(0.5,0.5)
  threshes[[13]] <- c(0.6,0.6)
  threshes
}