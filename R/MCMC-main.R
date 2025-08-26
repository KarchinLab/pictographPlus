#' main function to infer tumor evolution history
#' 
#' run MCMC chains to  infer the clonal evolution of tumors from single or multi-region sequencing data. 
#' This function automatically runs a pipeline of the tool. It models uncertainty of mutation cellular 
#' fraction (MCF) in small somatic mutations (SSMs) and copy number alterations (CNAs), assigning SSMs 
#' and CNAs to subclones using a Bayesian hierarchical model, and reconstruct tumor evolutionary trees 
#' that are constrained based on principles of lineage precedence, sum condition, and optionally by 
#' sample-presence. 
#' 
#' @param mutation_file a csv file that include information for SSMs. See vignette for details.
#' @param copy_number_file a csv file that include information for CNA. See vignette for details.
#' @param SNV_file a csv file that include information for germline heterozygous SNVs. See vignette for details.
#' @param outputDir output directory for saving all files.
#' @param sample_presence whether or not to use sample presence to separate the mutations; default: TRUE 
#' @param score scoring function to estimate the number of clusters. silhouette or BIC; default: silhuette
#' @param max_K user defined maximum number of clusters; default: 10
#' @param min_mutation_per_cluster minumum number of mutations in each cluster; default: 5
#' @param min_cluster_thresh minimum MCF for each cluster; default: 0.05
#' @param cluster_diff_thresh difference threshold to merge two clusters: default: 0.05
#' @param n.iter number of iterations by JAGS; default: 5000
#' @param n.burn number of burns by JAGS; default: 1000
#' @param thin number of thin by JAGS; default: 10
#' @param mc.cores number of cores to use for parallel computing; not applicable to windows; default: 8
#' @param inits additional parameters by JAGS.
#' @param LOH whether or not to include copy number segments that are copy neutral but LOH; default: FALSE
#' @param purity_min minimum purity for tumor samples; default: 0.2
#' @param driverFile list of driver genes used for visualization. See vignette for details.
#' @param cytobandFile list of cytoband regions used for visualization. See vignette for details.
#' @param alt_reads_thresh minimum number of alternative read count for a SSM to be included in the analysis; default: 0
#' @param vaf_thresh minimum VAF for a SSM to be included in the analysis; default: 0
#' @param tcn_normal_range range of total copy number considered as copy-neutral; default: c(1.75,2.3)
#' @param filter_cnv whether or not to filter copy number alterations; default: TRUE
#' @param smooth_cnv whether or not to process copy number alterations across samples to unify the segment start and end postions; default: TRUE
#' @param autosome to only include autosomes; default: TRUE
#' @param cnv_min_length minimum length of copy number alterations for it to be included in analysis
#' @param depth total_read counts to fill in if missing from the input
#' @export
runPictograph <- function(mutation_file,
                     copy_number_file=NULL,
                     SNV_file=NULL,
                     outputDir=NULL,
                     sample_presence=FALSE,
                     score="BIC", # either BIC or silhouette
                     max_K = 10, 
                     min_mutation_per_cluster=1, 
                     min_cluster_thresh=0.05,
                     cluster_diff_thresh=0.05,
                     n.iter=5000, 
                     n.burn=1000, 
                     thin=10, 
                     mc.cores=8, 
                     inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
                     LOH = FALSE,
                     purity_min=0.2,
                     driverFile = NULL,
                     cytobandFile = NULL,
                     alt_reads_thresh = 0, 
                     vaf_thresh = 0, 
                     cnv_min_length = 1000000,
                     depth=300,
                     tcn_normal_range = c(1.75, 2.3), 
                     filter_cnv = T, 
                     smooth_cnv = T, 
                     autosome = T,
                     dual_model=TRUE,
                     ploidy=2,
                     pval=0.05,
                     threshes=NULL
                     ) {
  
  data <- importFiles(mutation_file=mutation_file, 
                      copy_number_file=copy_number_file, 
                      outputDir=outputDir, 
                      SNV_file=SNV_file, 
                      LOH=LOH,
                      purity_min=purity_min,
                      alt_reads_thresh=alt_reads_thresh, 
                      vaf_thresh=vaf_thresh, 
                      cnv_min_length=cnv_min_length, 
                      tcn_normal_range=tcn_normal_range, 
                      filter_cnv=filter_cnv,
                      smooth_cnv=smooth_cnv,
                      autosome=autosome, 
                      pval=pval,
                      depth=depth)
  
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
  if (data$cnnull) { 
    
    if (sample_presence) {
      message("Mode: SSM only; using sample presence")
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
      sep_list <- separateMutationsBySamplePresence(input_data, min_mutation_per_cluster)
      
      # For each presence set, run clustering MCMC, calculate silhouette and BIC and choose best K
      all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                            min_cluster_thresh=min_cluster_thresh, cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                            n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
    } else {
      message("Mode: SSM only; not using sample presence")
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
                                            min_cluster_thresh=min_cluster_thresh, cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                            n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
    }
  } else {
    # using two step modeling
    if (dual_model) {
      if (sample_presence) {
        
        message("Mode: SSM and CNA; using sample presence")
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
        sep_list <- separateMutationsBySamplePresence(input_data, min_mutation_per_cluster)
        
        # For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
        all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              min_cluster_thresh=min_cluster_thresh, cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
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
        sep_list <- separateMutationsBySamplePresence(input_data, min_mutation_per_cluster)
        
        # For each presence set, run clustering MCMC
        all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=sample_presence, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              min_cluster_thresh=min_cluster_thresh, cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
        
      } else {
        
        message("Mode: SSM and CNA; not using sample presence")
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
        
        sep_list <- separateMutationsBySamplePresence(input_data, min_mutation_per_cluster)
        
        # For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
        all_set_results <- runMCMCForAllBoxes(sep_list, sample_presence=TRUE, ploidy=ploidy, max_K = max_K, min_mutation_per_cluster = min_mutation_per_cluster, 
                                              min_cluster_thresh=min_cluster_thresh, cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type1")
        # pick K
        set_k_choices <- writeSetKTable(all_set_results)
        
        # collect best chains
        # if (score == "silhouette") {
        #   best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$silhouette_K)
        # } else {
        #   best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$min_BIC)
        # }
        best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$min_BIC)
        
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
                                              min_cluster_thresh=min_cluster_thresh, cluster_diff_thresh = cluster_diff_thresh, inits = inits, 
                                              n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores, model_type = "type2")
        
       }
      
    } else {
      warning("Mode: SSM and CNA; using sample presence; single model NEED TO BE TESTED")
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
  if (sample_presence) {
    if (score=="silhouette") {
      best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$silhouette_K)
    } else {
      best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$BIC_K)
    }
  } else {
    if (score=="silhouette") {
      best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$silhouette_K)
    } else {
      best_set_chains <- collectBestKChains(all_set_results, chosen_K = set_k_choices$BIC_K)
    }
  }
  
  chains <- mergeSetChains(best_set_chains, input_data)
  chains <- assign("chains", chains, envir = .GlobalEnv)
  
  # plot MCMC tracing 
  png(paste(outputDir, "mcf.png", sep="/"))
  print(
    plotChainsMCF(chains$mcf_chain)
  )
  dev.off()
  
  # write mcf table
  mcfTable = writeClusterMCFsTable(chains$mcf_chain)
  colnames(mcfTable)=c("Cluster",c(colnames(data$y)))
  write.table(mcfTable, file=paste(outputDir, "mcf.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # write cluster assignment table
  clusterassignmentTable = writeClusterAssignmentsTable(chains$z_chain, Mut_ID = input_data$MutID)
  
  # record estimated icn and multiplicity information
  icnTable <- writeIcnTable(chains$icn_chain, Mut_ID = input_data$MutID)
  # write.table(icnTable, file=paste(outputDir, "icn_all.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  multiplicityTable <- writeMultiplicityTable(chains$m_chain, chains$icn_chain, Mut_ID = input_data$MutID)
  # write.table(multiplicityTable, file=paste(outputDir, "multiplicity_all.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  icnTableCN <- icnTable[data$is_cn==1,]
  multiplicityTableCN <- multiplicityTable[data$is_cn==1,]
  icnTableCN$Multiplicity <- multiplicityTableCN$Multiplicity
  write.table(icnTableCN, file=paste(outputDir, "CN_results.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # Clean copy number segments by removing segments with icn of 2 and multiplicity of 1
  toKeepIndex = c()
  for (i in seq_len(nrow(clusterassignmentTable))) {
    if (clusterassignmentTable[i,]$Mut_ID %in% icnTableCN$Mut_ID) {
      icnInfo <- icnTableCN %>% filter(Mut_ID==clusterassignmentTable[i,]$Mut_ID)
      if (icnInfo$icn==2 & icnInfo$Multiplicity==1) {
        toKeepIndex <- c(toKeepIndex, 0)
      } else {
        major_cn = max(icnInfo$Multiplicity, icnInfo$icn-icnInfo$Multiplicity)
        clusterassignmentTable[i,]$Mut_ID <- paste(clusterassignmentTable[i,]$Mut_ID, ";icn:", icnInfo$icn, ";", "major_cn:", major_cn, sep="")
        toKeepIndex <- c(toKeepIndex, 1)
      }
    } else {
      toKeepIndex <- c(toKeepIndex, 1)
    }
  }
  clusterassignmentTable$idx <- toKeepIndex
  clusterassignmentTable <- clusterassignmentTable %>% filter(toKeepIndex==1) %>% select(-idx)
  write.table(clusterassignmentTable, file=paste(outputDir, "clusterAssign.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  ###### generate combined cluster assignment output ###########################
  y <- as.data.frame(data$y) %>% rename_with(~ paste0(., "_variant_count"))
  n <- as.data.frame(data$n) %>% rename_with(~ paste0(., "_total_count"))
  tcn <- as.data.frame(data$tcn) %>% rename_with(~ paste0(., "_tcn"))
  tcn[tcn == 0] <- 2
  vaf <- data$y / data$n
  vaf <- as.data.frame(vaf) %>% rename_with(~ paste0(., "_vaf"))
  
  y <- rownames_to_column(y, var = "Mut_ID")
  n <- rownames_to_column(n, var = "Mut_ID")
  vaf <- rownames_to_column(vaf, var = "Mut_ID")
  tcn <- rownames_to_column(tcn, var = "Mut_ID")
  
  combined <- y %>%
    left_join(n, by = "Mut_ID") %>%
    left_join(tcn, by = "Mut_ID") %>%
    left_join(vaf, by = "Mut_ID") %>%
    select(order(names(.)))
  
  columns <- colnames(combined)
  mut_id_column <- "Mut_ID"
  other_columns <- setdiff(columns, mut_id_column)
  
  ordered_columns <- other_columns %>%
    data.frame(column_name = ., stringsAsFactors = FALSE) %>%
    mutate(
      sample = sub("_(variant_count|total_count|vaf|tcn)$", "", column_name), # Extract sample name
      suffix = sub("^.*_(variant_count|total_count|vaf|tcn)$", "\\1", column_name) # Extract suffix
    ) %>%
    mutate(
      suffix_order = case_when(
        suffix == "variant_count" ~ 1,
        suffix == "total_count"   ~ 2,
        suffix == "vaf"           ~ 3,
        suffix == "tcn"           ~ 4
      )
    ) %>%
    arrange(sample, suffix_order) %>%
    pull(column_name)
  
  combined <- combined %>%
    select(all_of(c(mut_id_column, ordered_columns)))
  
  clusterassignmentTable <- clusterassignmentTable %>%
    mutate(Mut_ID_processed = ifelse(
      grepl("^chr", Mut_ID),               # Check if Mut_ID starts with "chr"
      sub(";.*", "", Mut_ID),              # Extract the substring before the first semicolon
      Mut_ID                               # Otherwise, keep the original Mut_ID
    ))
  
  final_df <- combined %>%
    left_join(clusterassignmentTable, by = c("Mut_ID" = "Mut_ID_processed"))
  
  write.table(final_df, file=paste(outputDir, "mutationClusterAssign.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  #################################################################################
  
  # generate trees using different set of thresholds until at least one tree is available
  if (is.null(threshes)){
    threshes <- allThreshes()
  }
  for (thresh in threshes) {
    generateAllTrees(chains$mcf_chain, data$purity, lineage_precedence_thresh = thresh[1], sum_filter_thresh = thresh[2])
    if (length(all_spanning_trees) > 0) {
      break
    }
  }

  if (data$cnnull) {
    cncfTable <- data$cncf
  } else {
    cncfTable <- findCncf(data, input_data, chains)
  }
  scores <- calcTreeScores(chains$mcf_chain, all_spanning_trees, purity=data$purity)
  # scores <- calculateTreeScoreMutations(chains$mcf_chain, data, icnTable, cncfTable, multiplicityTable, clusterassignmentTable, data$purity, all_spanning_trees)
  scores <- assign("scores", scores, envir = .GlobalEnv)
  
  ##### find dirver mutations
  driverList = NULL
  filteredDriverList = NULL
  if (!is.null(driverFile)) {
    driverList <- read.csv(driverFile)
    filteredDriverList <- getDrivers(clusterassignmentTable, driverList, cytobandFile)
    write.table(filteredDriverList, file=paste(outputDir, "labelling.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  }
  
  # plot all possible trees
  plotAllTrees(outputDir, scores, all_spanning_trees, mcfTable, data, filteredDriverList)
  
  # highest scoring tree
  best_tree <- all_spanning_trees[[which(scores == max(scores))[length(which(scores == max(scores)))]]]
  write.table(best_tree, file=paste(outputDir, "tree.csv", sep="/"), quote = FALSE, sep = ",", row.names = F)
  
  # plot best and ensemble tree
  if (nrow(best_tree) >= 1) {
    png(paste(outputDir, "tree.png", sep="/"), width = 1200, height = 1600, res = 300)
    plotTree(best_tree, filteredDriverList, palette = viridis::viridis)
    dev.off()
  }
  
  # estimate purity
  cc <- best_tree %>% filter(parent=="root") %>% select(child)
  purity <- mcfTable %>% filter(Cluster %in% cc$child) %>% summarise(across(everything(), sum)) %>% select(-Cluster)
  # purity[purity>1] <- 1 # test
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

  # save all data
  save.image(file=paste(outputDir, "PICTographPlus.RData", sep="/"))
  
}

get_overlapping_gene <- function(chromosome, start1, end1, driverList) {
  overlap <- driverList %>%
    filter(chrom == chromosome & 
             start1 <= as.numeric(end) & 
             end1 >= as.numeric(start))
  if (nrow(overlap) > 0) {
    return(overlap$gene[1]) # Return the first overlapping gene
  } else {
    return(NA) # Return NA if no overlap is found
  }
}

getDrivers <- function(ClusterAssignmentTable, driverList, cytobandFile) {
  
  filtered_table <- ClusterAssignmentTable %>%
    mutate(
      # Check if Mut_ID starts with "chr"
      is_chr = str_starts(Mut_ID, "chr"),
      # Extract chrom, start, end, and icn from Mut_ID if it starts with "chr"
      chrom_start_end = ifelse(is_chr, str_extract(Mut_ID, "chr[0-9XY]+-[0-9]+-[0-9]+"), NA),
      icn_value = ifelse(is_chr, as.numeric(str_extract(Mut_ID, "(?<=icn:)\\d+")), NA)
    ) %>%
    rowwise() %>%
    filter(
      # Condition 1: Gene-based rows where the gene is in driverList
      (!is_chr & str_extract(Mut_ID, "^[^_]+") %in% driverList$gene) |
        # Condition 2: Segment-based rows where the segment contains a driver gene
        (is_chr & {
          # Extract chrom, start, and end from chrom_start_end
          chrom1 <- str_extract(chrom_start_end, "chr[0-9XY]+")
          start1 <- as.numeric(str_extract(chrom_start_end, "(?<=-)[0-9]+(?=-)"))
          end1 <- as.numeric(str_extract(chrom_start_end, "(?<=-)[0-9]+$"))
          
          # Check if any gene in driverList is within this segment
          matching_genes <- driverList %>%
            filter(chrom1 == chrom, start >= start1, end <= end1)
          
          # Further filter based on gene type and icn_value
          any(
            (matching_genes$gene_type == "oncogene" & icn_value > 2) |
              (matching_genes$gene_type == "tumor_suppressor" & icn_value < 2)
          )
        }) |
        # New Condition 3: If icn_value is 2, check if any driverList gene belongs to the segment and is in Mut_ID
        (is_chr & icn_value == 2 & {
          # Extract chrom, start, and end from chrom_start_end
          chrom1 <- str_extract(chrom_start_end, "chr[0-9XY]+")
          start1 <- as.numeric(str_extract(chrom_start_end, "(?<=-)[0-9]+(?=-)"))
          end1 <- as.numeric(str_extract(chrom_start_end, "(?<=-)[0-9]+$"))
          
          # Find matching genes from driverList within the segment
          matching_genes <- driverList %>%
            filter(chrom1 == chrom, start >= start1, end <= end1) %>%
            pull(gene)
          
          non_chr_genes <- ClusterAssignmentTable %>%
            filter(!str_starts(Mut_ID, "chr")) %>%
            pull(Mut_ID) %>%
            str_extract("^[^_]+")
          
          # Ensure at least one matching gene is mentioned in Mut_ID
          any(matching_genes %in% non_chr_genes)
        })
    )
  
  if (nrow(filtered_table) == 0) {
    return(NULL)
  }
  
  filtered_table <- filtered_table %>%
    rowwise() %>%
    mutate(
      # If Mut_ID_processed starts with "chr", process it
      Mut_ID = if (str_starts(Mut_ID_processed, "chr")) {
        # Extract chrom, start, and end from Mut_ID_processed
        chrom1 <- str_extract(Mut_ID_processed, "chr[0-9XY]+")
        start1 <- as.numeric(str_extract(Mut_ID_processed, "(?<=-)[0-9]+(?=-)"))
        end1 <- as.numeric(str_extract(Mut_ID_processed, "(?<=-)[0-9]+$"))
        
        # Find matching genes in driverList
        matching_genes <- driverList %>%
          filter(chrom == chrom1, start >= start1, end <= end1) %>%
          pull(gene)
        
        # Create the updated Mut_ID by appending genes and AMP/DEL
        gene_list <- if (length(matching_genes) > 0) paste(matching_genes, collapse = ";") else ""
        dup_del <- if (icn_value > 2) "AMP" else if (icn_value < 2) "DEL" else "LOH"
        paste(gene_list, dup_del, sep = "_")
        # paste(Mut_ID_processed, gene_list, dup_del, sep = "_")
      } else {
        # If Mut_ID_processed does not start with "chr", keep it as is
        Mut_ID_processed
      }
    ) %>%
    ungroup() %>%
    select(Mut_ID, Cluster, chrom_start_end)
  
  
  if (!is.null(cytobandFile)) {
    allowed_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
    cytobandTable = read.delim(cytobandFile, header = FALSE, stringsAsFactors = FALSE)
    cytobandTable <- cytobandTable[cytobandTable[[1]] %in% allowed_chromosomes, ]
    colnames(cytobandTable) <- c("chrom", "start", "end", "cytoband", "ext")
    cytobandTable <- unique(cytobandTable)
    
    filtered_table <- filtered_table %>%
      mutate(
        chromosome = sub("-.*", "", chrom_start_end),                    # Extract chromosome
        starts = as.numeric(sub("-.*", "", sub(".*?-", "", chrom_start_end))), # Extract starts
        ends = as.numeric(sub(".*-", "", chrom_start_end))                    # Extract ends
      )
    
    filtered_table <- filtered_table %>%
      rowwise() %>%
      mutate(
        cytoband = if (!is.na(chromosome)) {
          # Find the cytoband for the start position
          match_start <- which(cytobandTable$chrom == chromosome & 
                                 cytobandTable$start <= starts & 
                                 cytobandTable$end > starts)
          cytoband_start <- if (length(match_start) > 0) cytobandTable$cytoband[match_start[1]] else NA
          
          # Find the cytoband for the end position
          match_end <- which(cytobandTable$chrom == chromosome & 
                               cytobandTable$start <= ends & 
                               cytobandTable$end > ends)
          cytoband_end <- if (length(match_end) > 0) cytobandTable$cytoband[match_end[1]] else NA
          
          # Combine cytobands
          if (!is.na(cytoband_start) && cytoband_start == cytoband_end) {
            cytoband_start
          } else {
            paste(na.omit(c(cytoband_start, cytoband_end)), collapse = "-")
          }
        } else {
          NA  # If chromosome is NA, cytoband is NA
        },
        Mut_ID = if (!is.na(cytoband)) {
          paste(chromosome,cytoband, Mut_ID, sep = "_")
        } else {
          Mut_ID
        }
      ) %>%
      ungroup()
  } else {
    filtered_table <- filtered_table %>%
      rowwise() %>% mutate(
        Mut_ID = if (!is.na(chrom_start_end)) {
          paste(chrom_start_end, Mut_ID, sep = "_")
        } else {
          Mut_ID
        }
      ) %>%
      ungroup()
  }
  
  return(filtered_table)
}

#' Plot all trees with the highest scores
plotAllTrees <- function(outputDir, scores, all_spanning_trees, mcfTable, data, filteredDriverList) {
  # plot all tree with best scores
  
  outputDir = paste(outputDir, "all_trees", sep = "/")
  suppressWarnings(dir.create(outputDir))
  
  for (i in seq_len(length(which(scores == max(scores))))) {
    idx = which(scores == max(scores))[i]
    
    best_tree <- all_spanning_trees[[idx]]
    if (nrow(best_tree) >= 1 ) {
      write.table(best_tree, file=paste(outputDir, "/tree", i, ".csv", sep=""), quote = FALSE, sep = ",", row.names = F)
  
      png(paste(outputDir, "/tree", i, ".png", sep=""), width = 1600, height = 1200, res = 300)
      # plot tree
      plotTree(best_tree, filteredDriverList, palette = viridis::viridis)
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
  mTable <- writeMultiplicityTable(chains$m_chain, chains$icn_chain, Mut_ID = input_data$MutID)
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
  threshes[[2]] <- c(0.025,0.025)
  threshes[[3]] <- c(0.05,0.05)
  threshes[[4]] <- c(0.075,0.075)
  threshes[[5]] <- c(0.1,0.1)
  threshes[[6]] <- c(0.15,0.15)
  threshes[[7]] <- c(0.2,0.2)
  threshes[[8]] <- c(0.25,0.25)
  threshes[[9]] <- c(0.3,0.3)
  threshes[[10]] <- c(0.4,0.5)
  threshes[[11]] <- c(0.5,0.6)
  threshes[[12]] <- c(0.6,0.7)
  threshes[[13]] <- c(0.7,0.8)
  threshes[[14]] <- c(0.8,1)
  threshes
}
