#' read input data file and store in PICTOGRAPH input format
#' @export
#' @param mutation_file mutation data file that contains columns "sample", "mutation", "chrom", "start", "end", total_reads", and "alt_reads";
#' @param copy_number_file copy number file that contains columns "sample", "chrom", "start", "end", "tcn"
importFiles <- function(mutation_file, 
                        copy_number_file=NULL, 
                        SNV_file=NULL, 
                        outputDir=NULL,
                        stat_file=NULL, 
                        alt_reads_thresh = 0, 
                        vaf_thresh = 0, 
                        cnv_max_dist=2000, 
                        cnv_max_percent=0.30, 
                        tcn_normal_range=c(1.8, 2.2), 
                        smooth_cnv=F, 
                        autosome=T, 
                        mc.cores=8, 
                        pval=0.05) {
  
  if (is.null(outputDir)) {
    outputDir = getwd()
  }
  
  if (!is.null(copy_number_file)) {
    # keep mutations if alt_reads >= alt_reads_thresh and vaf >= vaf_thresh
    mutation_data = importMutationFile(mutation_file, alt_reads_thresh, vaf_thresh)
    name_order <- colnames(mutation_data$y)
    
    copy_number_data = importCopyNumberFile(copy_number_file, 
                                            outputDir, 
                                            SNV_file, 
                                            stat_file, 
                                            name_order,
                                            cnv_max_dist, 
                                            cnv_max_percent, 
                                            tcn_normal_range, 
                                            smooth_cnv, 
                                            autosome, 
                                            pval, 
                                            mc.cores)
    
    mutation_data$tcn = copy_number_data$tcn[, name_order, drop=FALSE]
    
    # bind SSM and CNA
    mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)), rep(1,nrow(copy_number_data$tcn)))
    mutation_data$y <- rbind(mutation_data$y, copy_number_data$tcn_alt[, name_order, drop=FALSE])
    mutation_data$n <- rbind(mutation_data$n, copy_number_data$tcn_tot[, name_order, drop=FALSE])
    mutation_data$MutID <- c(mutation_data$MutID, rownames(copy_number_data$tcn))
    mutation_data$I <- mutation_data$I + nrow(copy_number_data$tcn)
    
    mutation_data$overlap = resolveOverlap(mutation_data)
    
    mutation_data$tcn <- mutation_data$overlap %*% mutation_data$tcn
    
    
    overlap <- mutation_data$overlap
    
    colnames(overlap) <- sapply(colnames(mutation_data$overlap), function(col) {which(rownames(mutation_data$overlap)==col)})
    
    q <- vector("numeric", nrow(overlap))
    for (i in 1:nrow(overlap)) {
      q[i] <- ifelse(length(which(overlap[i,] == 1)) > 0, as.numeric(names(which(overlap[i,] == 1))[1]),i)
    }
    mutation_data$q <- q
    
  } else {
  
    mutation_data = importMutationFileOnly(mutation_file, alt_reads_thresh, vaf_thresh)
    # name_order <- colnames(mutation_data$y)
    
    mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)))
  }
  
  
  
  # warning("resolveOverlap: need to check whether a mutation overlaps with two CNA segs")

  return(mutation_data)
}

#' check whether a mutation overlaps a CNA region
resolveOverlap <- function(mutation_data) {
  mut_count = nrow(mutation_data$y)
  cna_count = nrow(mutation_data$tcn)
  output = matrix(0, nrow=mut_count, ncol=cna_count)
  rownames(output) = rownames(mutation_data$y)
  colnames(output) = rownames(mutation_data$tcn)
  # print(output)
  cnas = rownames(mutation_data$tcn)
  mutations = mutation_data$position
  for (i in seq_len(mut_count)) {
    for (j in seq_len(cna_count)) {
      if (mutation_data$is_cn[i]==0) {
        chrSplit = strsplit(cnas[j], split='-')[[1]]
        chr = chrSplit[1]
        start_pos = strtoi(chrSplit[2])
        end_pos = strtoi(chrSplit[3])
        
        if (chr==mutations[i,2]) {
          if (strtoi(mutations[i,3]) >= start_pos && strtoi(mutations[i,4]) <= end_pos) {
            output[i,j] = 1
            break
          }
        }
      } else {
        if (rownames(mutation_data$y)[i]==rownames(mutation_data$tcn)[j]) {
          output[i,j] = 1
        }
      }
    }
  }
  return(output)
}

#' merge segments with similar coordinates
#' the maximum allowed distance between the start or end position of two segments is the max(cnv_max_dist, min(length of two segments)*cnv_max_percent)
smoothCNV <- function(data, cnv_max_dist=2000, cnv_max_percent=0.30) {
  # smooth_data <- as_tibble(data.frame(matrix(nrow=0, ncol=ncol(data))))
  # colnames(smooth_data) <- colnames(data)
  
  # samples <- unique(data$sample)
  
  message("improvements needed for smoothCNV under preprocessing.R; currently detecting two segments as the same if overlap and star/end position within certain distance; merge the two segments as the outcome")
  data <- data %>% 
    mutate(numeric_chrom = as.numeric(sub("chr", "", chrom))) %>% 
    arrange(numeric_chrom, start) %>%
    select(-numeric_chrom)
  
  indexList = seq_len(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    if (i %in% indexList) {
      # print(data[i,])
      max_dis = max(cnv_max_dist, cnv_max_percent*(data[i,]$end-data[i,]$start))
      index_selected = which((data$chrom==data[i,]$chrom)&(data$start<=data[i,]$end)&(data[i,]$start<=data$end)&(abs(data$start-data[i,]$start)<max_dis)&(abs(data$end-data[i,]$end)<max_dis))
      # index_selected = which((data$chrom==data[i,]$chrom)&(abs(data$start-data[i,]$start)<max_dis)&(abs(data$end-data[i,]$end)<max_dis))
      if (length(index_selected) > 1) {
        # list of cnv segments able to merge
        cnv_selected = data[index_selected,]
        # start = min(cnv_selected$start)
        # end = max(cnv_selected$end)
        data[index_selected,]$start = min(cnv_selected$start)
        data[index_selected,]$end = max(cnv_selected$end)
        data[index_selected,]$drivers = paste(unique(unlist(strsplit(cnv_selected$drivers, ";"))), collapse=";")
        data[index_selected,]$genes = paste(unique(unlist(strsplit(cnv_selected$genes, ";"))), collapse=";")
      }
      
      indexList <- indexList[!(indexList %in% index_selected)]
    }
  }
  
  return(data)
}

#' import copy number file
#' @param cnv_max_dist: maximum of distance allowed between two segments to assign as the same one
#' @param cnv_max_percent: maximum percentage of distance allowed between two segments to assign as the same one
#' @param smooth_cnv: process input CNV to merge  segments with similar distance
importCopyNumberFile <- function(copy_number_file, outputDir, SNV_file=NULL, stat_file=NULL, name_order=NULL, cnv_max_dist=2000, cnv_max_percent=0.30, tcn_normal_range=c(1.8, 2.2), smooth_cnv=F, autosome=T, pval=0.05, mc.cores=8) {
  
  data <- read_csv(copy_number_file, show_col_types = FALSE) # read copy number csv file
  
  if ("baf" %in% colnames(data)) {
    message("inferring allele-specific copy number using BAF")
  } else if (!is.null(SNV_file)) {
    message("inferring allele-specific copy number using heterozygous SNVs")
    data <- check_sample_LOH(data, outputDir, SNV_file, tcn_normal_range=tcn_normal_range, pval=pval) # check unimodality at both normal and tumor sample
    data <- data[data$to_keep==1,] # keep rows is to_keep is 1
  } else {
    stop("Please provide either a baf column in the copy number file or supply a SNV file with heterozygous SNV counts")
  }
  
  
  data$tcn[data$tcn==2] <- 2.01 # make tcn=2 -> 2.01 to avoid confusion during sample presence
  
  # Smooth segments so segments with different start/end position are treated the same, 
  # conditioned on overlapping and distance between positions
  if (smooth_cnv) {
    data <- smoothCNV(data, cnv_max_dist, cnv_max_percent)
  }
  
  data$CNA = paste(data$chrom, data$start, data$end, sep = '-')
  
  # keep only autosome
  if (autosome) {
    data <- data[!(data$chrom %in% c('chrX', 'chrY', 'chrM')),]
  }
  
  # convert to cna by sample matrix 
  output_data <- list()
  output_data = as.matrix(data[c("sample", "CNA", "tcn")] %>% pivot_wider(names_from = sample, values_from = tcn, values_fill = 2))
  
  # fill missing CNA with 2
  message("FILL MISSING TCN WITH 2; NEED TO FIX FOR CHRX/Y FOR VALUES_FILL in importCopyNumberFile in preprocessing.R")
  rownames(output_data) <- output_data[,'CNA']
  output_data <- output_data[,-1, drop=FALSE]
  rowname = rownames(output_data)
  colname = colnames(output_data)
  output_data <- matrix(as.numeric(output_data), ncol = ncol(output_data))
  rownames(output_data) = rowname
  colnames(output_data) = colname
  
  output_data <- add_missing_column(name_order, output_data, 2)
  
  
  if ("baf" %in% colnames(data)) {
    baf = as.matrix(data[c("sample", "CNA", "baf")] %>% pivot_wider(names_from = sample, values_from = baf, values_fill = 0.5))
    rownames(baf) <- baf[,'CNA']
    baf <- baf[,-1, drop=FALSE]
    rowname = rownames(baf)
    colname = colnames(baf)
    baf <- matrix(as.numeric(baf), ncol = ncol(baf))
    rownames(baf) = rowname
    colnames(baf) = colname
    
    tcn_tot <- matrix(200, nrow(output_data), ncol(output_data))
    rownames(tcn_tot) <- rownames(output_data)
    colnames(tcn_tot) <- colnames(output_data)
    
    tcn_tot <- add_missing_column(name_order, tcn_tot, 200)
    
    tcn_alt <- matrix(round(tcn_tot * baf), nrow(output_data),ncol(output_data))
    rownames(tcn_alt) <- rownames(output_data)
    colnames(tcn_alt) <- colnames(output_data)
    
    tcn_alt <- add_missing_column(name_order, tcn_alt, 100)
    
  } else if (!is.null(SNV_file)) {
    tcn_ref <- list()
    tcn_ref = as.matrix(data[c("sample", "CNA", "tcn_ref")] %>% pivot_wider(names_from = sample, values_from = tcn_ref, values_fill = 0))
    rownames(tcn_ref) <- tcn_ref[,'CNA']
    tcn_ref <-tcn_ref[,-1, drop=FALSE]
    rowname = rownames(tcn_ref)
    colname = colnames(tcn_ref)
    tcn_ref <- matrix(as.numeric(tcn_ref), ncol = ncol(tcn_ref))
    rownames(tcn_ref) = rowname
    colnames(tcn_ref) = colname
    
    tcn_alt <- list()
    tcn_alt = as.matrix(data[c("sample", "CNA", "tcn_alt")] %>% pivot_wider(names_from = sample, values_from = tcn_alt, values_fill = 0))
    rownames(tcn_alt) <- tcn_alt[,'CNA']
    tcn_alt <-tcn_alt[,-1, drop=FALSE]
    rowname = rownames(tcn_alt)
    colname = colnames(tcn_alt)
    tcn_alt <- matrix(as.numeric(tcn_alt), ncol = ncol(tcn_alt))
    rownames(tcn_alt) = rowname
    colnames(tcn_alt) = colname
    
    tcn_tot <- tcn_ref + tcn_alt
    alt_mean <- round(mean(tcn_tot)/2)
    tot_mean <- round(mean(tcn_tot))
    tcn_alt[tcn_tot==0] <- alt_mean
    tcn_tot[tcn_tot==0] <- tot_mean
    
    tcn_tot <- add_missing_column(name_order, tcn_tot, tot_mean)
    tcn_alt <- add_missing_column(name_order, tcn_alt, alt_mean)
  }
  
  
  # # check lOH using HET SNVs; 
  # # if both germline and tumor SNV counts provided, will check LOH
  # # otherwise, keep CNA outside normal range only
  # if (!is.null(SNV_file)) {
  #   # if there is a CN_stat file already, load
  #   if (!is.null(stat_file)) {
  #     temp_data <- read_csv(stat_file)
  #     samples <- colnames(output_data)
  #     to_keep_index <- suppressWarnings(which(as.numeric(temp_data$tumor_pval)<pval & as.numeric(temp_data$germline_pval)>pval))
  #   } else {
  #     to_keep_index = check_LOH(output_data, outputDir, SNV_file, tcn_normal_range, pval, mc.cores, sim_iter)
  #   }
  # } else {
  #   # to_keep_index = which(rowSums(output_data<tcn_normal_range[1]|output_data>tcn_normal_range[2])>0)
  #   errorCondition("Current implementation requires a germline heterozygous SNV file, please provide one; Will be updating w/o SNV file")
  #   stop()
  # }
  
  # tcn_ref <- read_csv(paste(outputDir, "tcf_ref.csv", sep="/"))
  # tcn_ref_matrix <- as.matrix(tcn_ref[,!names(tcn_ref) %in% "chrom", drop = FALSE])
  # rownames(tcn_ref_matrix) <- tcn_ref[["chrom"]]
  # tcn_ref = tcn_ref_matrix[to_keep_index,, drop=FALSE]
  # 
  # tcn_alt <- read_csv(paste(outputDir, "tcf_alt.csv", sep="/"))
  # tcn_alt_matrix <- as.matrix(tcn_alt[,!names(tcn_alt) %in% "chrom", drop = FALSE])
  # rownames(tcn_alt_matrix) <- tcn_alt[["chrom"]]
  # tcn_alt = tcn_alt_matrix[to_keep_index,, drop=FALSE]
  # 
  # output_data = output_data[to_keep_index,,drop=FALSE]
  # output_data[output_data==2] <- 2.01
  
  return_data <- list()
  return_data$tcn <- output_data
  return_data$tcn_tot <- tcn_tot
  return_data$tcn_alt <- tcn_alt
  return(return_data)
}

add_missing_column <- function(name_order, output_data, val) {
  matches <- name_order %in% colnames(output_data)
  sorted_matrix <- output_data[, colnames(output_data) %in% name_order, drop = FALSE]
  sorted_matrix <- sorted_matrix[, match(name_order[matches], colnames(output_data)), drop = FALSE]
  non_matches <- name_order[!matches]
  if (length(non_matches) > 0) {
    zero_cols <- matrix(val, nrow = nrow(output_data), ncol = length(non_matches))
    colnames(zero_cols) <- non_matches
    sorted_matrix <- cbind(sorted_matrix, zero_cols)
  }
  output_data <- sorted_matrix
  
  return(output_data)
}

#' @import LaplacesDemon
#' @import parallel
#' @import diptest
check_sample_LOH <- function(data, outputDir, SNV_file, tcn_normal_range=c(1.52, 2.46), pval=0.05) {
  # print(output_data)
  
  # tcn_normal_range=c(1.8, 2.2)
  # pval=0.05
  SNV_data <- read_csv(SNV_file)
  samples <- unique(data$sample)
  # output <- rownames_to_column(as.data.frame(output_data), var="chrom") 
  # output <- output %>% separate(col="chrom", into = c("chrom", "start", "end"), sep="-")
  # output <- output %>% mutate(start = as.integer(start), end = as.integer(end))
  # 
  pdf(paste(outputDir, "sample_modality.pdf", sep="/"), width = 12, height = 18)
  par(mfrow=c(6,3), cex.main = 1, cex.axis = 1, cex.lab = 1, cex = 1)
  
  # rowLen =10
  rowLen = nrow(data)
  
  # output data from this function
  tcn_ref <- rep(0, nrow(data))
  tcn_alt <- rep(0, nrow(data))
  
  # output_stats <- matrix(, nrow=0, ncol=8+length(samples))
  # colnames(output_stats) = c("chrom", "germline_depth", "germline_unimodal_prob", "germline_unimodality", "germline_pval", "tumor_depth", "tumor_unimodal_prob", samples, "tumor_pval")
  
  to_keep_index = c()
  # pval_list = c()
  # set.seed(123)
  for (i in 1:rowLen) {
    # i = 1
    # print(i)
    toSkip = FALSE
    SNV_temp <- SNV_data %>% filter(chroms==data[i,]$chrom & position>=data[i,]$start & position<=data[i,]$end)
    sample <- data[i,]$sample
    # SNV_temp
    if (nrow(SNV_temp) > 2) {
      # for a given segment, if tcn in all samples within tcn_normal_range, assuming copy neutral and test of LOH
      # using unimodal
      # vaf_list <- list()
      # um_list <- list()
      # 
      # germline_depth = round(sum(SNV_temp$germline_alt + SNV_temp$germline_ref) / nrow(SNV_temp))
      # results <- mclapply(1:sim_iter, function(i) {simulation_hidden(mcf=0, icn=2, minor_cn=1, depth=germline_depth, num_SNV=nrow(SNV_temp), seed=i, normal_prop=0)}, mc.cores=mc.cores)
      # um_prob <- sum(unlist(results)==TRUE) / sim_iter
      # 
      vaf_germline <- SNV_temp$germline_alt / (SNV_temp$germline_alt + SNV_temp$germline_ref)
      germline_test <- dip.test(vaf_germline)
      # vaf_list[[length(vaf_list) + 1]] <- vaf_germline
      # um_list[[length(um_list) + 1]] <- is.unimodal(vaf_germline)
      # 
      # tumor_depth = round(sum(SNV_temp[, 7:ncol(SNV_temp)])/(nrow(SNV_temp)*length(samples)))
      # results <- mclapply(1:sim_iter, function(i) {simulation_hidden(mcf=0, icn=2, minor_cn=1, depth=tumor_depth, num_SNV=nrow(SNV_temp), seed=i, normal_prop=0)}, mc.cores=mc.cores)
      # um_prob_tumor <- sum(unlist(results)==TRUE) / sim_iter

      alt = paste(sample, "alt", sep="_")
      ref = paste(sample, "ref", sep="_")
      
      vaf <- SNV_temp[[alt]] / (SNV_temp[[alt]] + SNV_temp[[ref]])

      # vaf_list[[length(vaf_list) + 1]] <- vaf
      # um_list[[length(um_list) + 1]] <- is.unimodal(vaf)
      
      tumor_test <- dip.test(vaf)
      
      if (is.unimodal(vaf)) { # if similar to germline distribution, alt and ref will be the sum of all
        tcn_alt[i] = round(mean(SNV_temp[[alt]]))
        tcn_ref[i] = round(mean(SNV_temp[[ref]]))
      } else if (is.trimodal(vaf)) { # else, cluster the vaf into two and take one cluster
        kmeans_result <- kmeans(vaf, centers = 3)
        cluster_number = which.max(kmeans_result$centers)
        tcn_alt[i] = round(mean(SNV_temp[[alt]][kmeans_result$cluster == cluster_number]))
        tcn_ref[i] = round(mean(SNV_temp[[ref]][kmeans_result$cluster == cluster_number]))
      } else if (is.bimodal(vaf)) {
        kmeans_result <- kmeans(vaf, centers = 2)
        cluster_number = which.max(kmeans_result$centers)
        tcn_alt[i] = round(mean(SNV_temp[[alt]][kmeans_result$cluster == cluster_number]))
        tcn_ref[i] = round(mean(SNV_temp[[ref]][kmeans_result$cluster == cluster_number]))
      } else {
        tcn_alt[i] = -1
        tcn_ref[i] = -1
        toSkip = TRUE
      }

      # output_stats <- rbind(output_stats, matrix(c(rownames(output_data)[i], germline_depth, um_prob,um_list[[1]], germline_test$p.value, tumor_depth, um_prob_tumor, unlist(um_list[2:length(um_list)]), tumor_test$p.value), nrow=1))
      # pval_list <- c(pval_list, tumor_test$p.value)
      # check if a segment is within tcn_normal_range for LOH
      if (toSkip) {
        to_keep_index <- c(to_keep_index, 0)
      } else {
        if (!germline_test$p.value < pval/nrow(data)) {
          if (data[i,]$tcn > tcn_normal_range[1] & data[i,]$tcn < tcn_normal_range[2]) {
            if (tumor_test$p.value < pval/nrow(data)) {
              
              # to_keep_index <- c(to_keep_index, 1)
              warning("check_sample_LOH diptest seems not working; not detecting LOH; more testing needed")
              
              to_keep_index <- c(to_keep_index, 0)
              plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, ", pval: ", tumor_test$p.value, sep=""))
              # print(i)
              # stop()
            } else {
              to_keep_index <- c(to_keep_index, 0)
            }
          } else {
            to_keep_index <- c(to_keep_index, 1)
            plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, sep=""))
          }
        } else {
          to_keep_index <- c(to_keep_index, 0)
        }
      }
      
    } else {
      to_keep_index <- c(to_keep_index, 0)
      # message("what if no SNV in the region; DROP")
      # output_stats <- rbind(output_stats, matrix(c(rownames(output_data)[i], rep("N/A", 7+length(samples))), nrow=1))
    }
  }
  dev.off()
  
  data <- data %>% add_column(to_keep = to_keep_index, tcn_ref = tcn_ref, tcn_alt = tcn_alt)
  
  return(data)
}

#' @import LaplacesDemon
#' @import parallel
check_LOH <- function(output_data, outputDir, SNV_file, tcn_normal_range=c(1.8, 2.2), pval=0.05, mc.cores=8, sim_iter=100) {
  # print(output_data)
  
  # tcn_normal_range=c(1.8, 2.2)
  # pval=0.05
  SNV_data <- read_csv(SNV_file)
  samples <- colnames(output_data)
  output <- rownames_to_column(as.data.frame(output_data), var="chrom") 
  output <- output %>% separate(col="chrom", into = c("chrom", "start", "end"), sep="-")
  output <- output %>% mutate(start = as.integer(start), end = as.integer(end))
  
  pdf(paste(outputDir, "modality.pdf", sep="/"), width = (length(samples)+1)*6, height = (length(samples)+1)*9)
  par(mfrow=c(length(samples)*3,(length(samples)+1)*2), cex.main = 1, cex.axis = 1, cex.lab = 1, cex = 1)
  
  # rowLen =10
  rowLen = nrow(output_data)
  
  # output data from this function
  tcn_ref = output_data
  tcn_ref[,] <- 0
  tcn_alt = output_data
  tcn_alt[,] <- 0
  
  output_stats <- matrix(, nrow=0, ncol=8+length(samples))
  colnames(output_stats) = c("chrom", "germline_depth", "germline_unimodal_prob", "germline_unimodality", "germline_pval", "tumor_depth", "tumor_unimodal_prob", samples, "tumor_pval")
  
  to_keep_index = c()
  
  # set.seed(123)
  for (i in 1:rowLen) {
    # 326; 316
    # set.seed(123)

    # i = 326
    # print(output[i,])
    SNV_temp <- SNV_data %>% filter(chroms==output[i,]$chrom & position>=output[i,]$start & position<=output[i,]$end)
    # SNV_temp
    if (nrow(SNV_temp) > 2) {
      # for a given segment, if tcn in all samples within tcn_normal_range, assuming copy neutral and test of LOH
      # using unimodal
      vaf_list <- list()
      um_list <- list()
      
      germline_depth = round(sum(SNV_temp$germline_alt + SNV_temp$germline_ref) / nrow(SNV_temp))
      # germline_depth
      # nrow(SNV_temp)
      
      # Estimate the probability of multi-modality given num SNVs and depth
      # start_time <- proc.time()
      # sim_iter = 100 # number of simulations
      
      # sim_unimodal = 0
      # for (iter in seq_len(sim_iter)) {
      #   sim_one = simulation_hidden(mcf=0, icn=2, minor_cn=1, depth=germline_depth, num_SNV=nrow(SNV_temp), seed=iter, normal_prop=0)
      #   if (sim_one) {
      #     sim_unimodal = sim_unimodal + 1
      #   }
      # }
      # print(sim_unimodal)
      # um_prob <- sim_unimodal/sim_iter
      # mm_prob <- 1 - sim_unimodal/sim_iter
      # mm_prob
      results <- mclapply(1:sim_iter, function(i) {simulation_hidden(mcf=0, icn=2, minor_cn=1, depth=germline_depth, num_SNV=nrow(SNV_temp), seed=i, normal_prop=0)}, mc.cores=mc.cores)
      um_prob <- sum(unlist(results)==TRUE) / sim_iter
      # mm_prob <- 1 - um_prob
      
      # end_time <- proc.time()
      # time_taken <- end_time - start_time
      # print(time_taken)
      
      vaf_germline <- SNV_temp$germline_alt / (SNV_temp$germline_alt + SNV_temp$germline_ref)
      vaf_list[[length(vaf_list) + 1]] <- vaf_germline
      um_list[[length(um_list) + 1]] <- is.unimodal(vaf_germline)
      # the probability of seeing (non)unimodal in germline SNV
      # if (ifelse(is.unimodal(vaf_germline), 1-mm_prob, mm_prob) >= pval) {
      # } 
      
      # vaf_germline <- SNV_temp$germline_alt / (SNV_temp$germline_alt + SNV_temp$germline_ref)
      # title = paste(output[i,]$chrom, output[i,]$start, output[i,]$end, sep="_")
      # plot(density(vaf_germline), xlim=c(0,1), main = title)
      # qqnorm(vaf_germline, main = title)
      # qqline(vaf_germline, col="grey")
      
      # ks_list = c()
      tumor_depth = round(sum(SNV_temp[, 7:ncol(SNV_temp)])/(nrow(SNV_temp)*length(samples)))
      results <- mclapply(1:sim_iter, function(i) {simulation_hidden(mcf=0, icn=2, minor_cn=1, depth=tumor_depth, num_SNV=nrow(SNV_temp), seed=i, normal_prop=0)}, mc.cores=mc.cores)
      um_prob_tumor <- sum(unlist(results)==TRUE) / sim_iter
      # mm_prob_tumor <- 1 - um_prob_tumor
      
      for (j in seq_len(length(samples))) {
        # j = 1
        
        alt = paste(samples[j], "alt", sep="_")
        ref = paste(samples[j], "ref", sep="_")
        
        # tumor_depth = round(sum(SNV_temp[[alt]] + SNV_temp[[ref]]) / nrow(SNV_temp))
        
        vaf <- SNV_temp[[alt]] / (SNV_temp[[alt]] + SNV_temp[[ref]])
        # title = paste(output[i,]$chrom, output[i,]$start, output[i,]$end, samples[j], sep="_")
        # title = samples[j]
        
        # ks = ks.test(vaf_germline, vaf)
        # ks = is.unimodal(vaf)
        vaf_list[[length(vaf_list) + 1]] <- vaf
        um_list[[length(um_list) + 1]] <- is.unimodal(vaf)
        
        if (is.unimodal(vaf)) { # if similar to germline distribution, alt and ref will be the sum of all
          tcn_alt[i,j] = round(mean(SNV_temp[[alt]]))
          tcn_ref[i,j] = round(mean(SNV_temp[[ref]]))
          # plot(density(vaf), xlim=c(0,1), main = title)
          # qqnorm(vaf, main = title)
          # qqline(vaf, col="green")
        } else if (is.trimodal(vaf)) { # else, cluster the vaf into two and take one cluster
          kmeans_result <- kmeans(vaf, centers = 3)
          cluster_number = which.max(kmeans_result$centers)
          tcn_alt[i,j] = round(mean(SNV_temp[[alt]][kmeans_result$cluster == cluster_number]))
          tcn_ref[i,j] = round(mean(SNV_temp[[ref]][kmeans_result$cluster == cluster_number]))
          # plot(density(vaf), xlim=c(0,1), main = title)
          # qqnorm(vaf, main = title)
          # qqline(vaf, col="red")
        } else {
          kmeans_result <- kmeans(vaf, centers = 2)
          cluster_number = which.max(kmeans_result$centers)
          tcn_alt[i,j] = round(mean(SNV_temp[[alt]][kmeans_result$cluster == cluster_number]))
          tcn_ref[i,j] = round(mean(SNV_temp[[ref]][kmeans_result$cluster == cluster_number]))
        }
        # ks_list <- c(ks_list, ks)
      }
      
      germline_test <- binom.test(x=sum(unlist(um_list)[1]==TRUE), n=1, p=um_prob)
      # germline_test$p.value # if smaller then pval, disgard the segment (probably gemline CNA event)
      
      tumor_test <- binom.test(x=sum(unlist(um_list)[2:length(um_list)]==TRUE), n=length(um_list)-1, p=um_prob_tumor)
      # tumor_test$p.value # if smaller then pval, keep the CNA (significantCNA event)
      # print(ks_list)
      
      # um_list[[length(um_list) + 1]] <- germline_depth
      # um_list[[length(um_list) + 1]] <- um_prob
      # um_list[[length(um_list) + 1]] <- germline_test$p.value
      # um_list[[length(um_list) + 1]] <- tumor_test$p.value
      # coln = rownames(output_data)[i]
      # output_stats <- list(unlist(um_list))
      output_stats <- rbind(output_stats, matrix(c(rownames(output_data)[i], germline_depth, um_prob,um_list[[1]], germline_test$p.value, tumor_depth, um_prob_tumor, unlist(um_list[2:length(um_list)]), tumor_test$p.value), nrow=1))
      # check if a segment is within tcn_normal_range for LOH
      # if (all(output[i,4:ncol(output)] > tcn_normal_range[1] & output[i,4:ncol(output)] < tcn_normal_range[2])) {
      #   if (all(!ks_list)) {
      #     to_keep_index <- c(to_keep_index, i)
      #   } 
      # } else {
      #   to_keep_index <- c(to_keep_index, i)
      # }
      
      # if (germline_test$p.value > pval & tumor_test$p.value < pval) {
      if (is.unimodal(vaf_germline) & tumor_test$p.value < pval) {
        to_keep_index <- c(to_keep_index, i)
        for (idx in seq_len(length(vaf_list))) {
          vaf <- vaf_list[[idx]]
          color = "red"
          if (um_list[[idx]]) {color = "green"}
          if (idx == 1) {
            title = paste(output[i,]$chrom, output[i,]$start, output[i,]$end, sep="_")
          } else {
            title = samples[idx-1]
          }
          # print(title)
          plot(density(vaf), xlim=c(0,1), main = title)
          qqnorm(vaf, main = ifelse(idx==1, paste("pval: ", tumor_test$p.value, sep=""), paste("tcn: ", output[i, idx+2], sep="")))
          qqline(vaf, col=color)
        }
      }
    } else {
      # message("what if no SNV in the region; DROP")
      output_stats <- rbind(output_stats, matrix(c(rownames(output_data)[i], rep("N/A", 7+length(samples))), nrow=1))
    }
  }
  dev.off()
  
  output_stats <- as.data.frame(output_stats)
  tcn_ref <- as.data.frame(tcn_ref)
  tcn_ref$chrom <- row.names(tcn_ref)
  tcn_alt <- as.data.frame(tcn_alt)
  tcn_alt$chrom <- row.names(tcn_alt)
  # output_stats
  # print(length(to_keep_index))
  write_csv(output_stats, file=paste(outputDir, "CN_stats.csv", sep="/"))
  write_csv(tcn_ref, file=paste(outputDir, "tcf_ref.csv", sep="/"))
  write_csv(tcn_alt, file=paste(outputDir, "tcf_alt.csv", sep="/"))
  return(to_keep_index)
}

#' import mutation file
importMutationFile <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0) {
  data <- read_csv(mutation_file, show_col_types = FALSE)
  data <- data %>% filter(alt_reads/total_reads>=vaf_thresh, alt_reads>=alt_reads_thresh)
  output_data <- list()

  output_data$y <- as.matrix(data[c("mutation", "sample", "alt_reads")] %>% pivot_wider(names_from = sample, values_from = alt_reads, values_fill = 0))
  rownames(output_data$y) <- output_data$y[,'mutation']
  output_data$y <- output_data$y[,-1, drop=FALSE]
  rowname = rownames(output_data$y)
  colname = colnames(output_data$y)
  output_data$y <- matrix(as.numeric(output_data$y), ncol = ncol(output_data$y))
  rownames(output_data$y) = rowname
  colnames(output_data$y) = colname

  output_data$MutID <- rowname

  output_data$n <- as.matrix(data[c("mutation", "sample", "total_reads")] %>% pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0))
  rownames(output_data$n) <- output_data$n[,'mutation']
  output_data$n <- output_data$n[,-1, drop=FALSE]
  rowname = rownames(output_data$n)
  colname = colnames(output_data$n)
  output_data$n <- matrix(as.numeric(output_data$n), ncol = ncol(output_data$n))
  rownames(output_data$n) = rowname
  colnames(output_data$n) = colname

  if ("purity" %in% colnames(data)) {
    output_data$purity <- as.matrix(data[c("mutation", "sample", "purity")] %>% pivot_wider(names_from = sample, values_from = purity, values_fill = 0))
    rownames(output_data$purity) <- output_data$purity[,'mutation']
    output_data$purity <- output_data$purity[,-1, drop=FALSE]
    rowname = rownames(output_data$purity)
    colname = colnames(output_data$purity)
    output_data$purity <- matrix(as.numeric(output_data$purity), ncol = ncol(output_data$purity))
    rownames(output_data$purity) = rowname
    colnames(output_data$purity) = colname
    output_data$purity = colSums(output_data$purity) / colSums(!!output_data$purity)
  } else {
    output_data$purity = rep(0.8, ncol(output_data$y))
  }
  
  if (any((output_data$y - output_data$n) > 0)) {
    warning("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
    stop()
  }

  if (any(output_data$n==0)) {
    print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- round(mean(output_data$n))
  }

  output_data$S = ncol(output_data$y)
  output_data$I = nrow(output_data$y)

  mutation_position = unique(data[c("mutation", "chrom", "start", "end")])
  if (!all(sort(unique(data[c("mutation", "chrom", "start", "end")]) %>% pull(mutation)) == sort(output_data$MutID))) {
    warning("Some mutations may have duplicated chromosome information; keeping the first occurence.")
    mutation_position = mutation_position[match(unique(mutation_position$mutation), mutation_position$mutation), ]
  }

  output_data$position = mutation_position
  return(output_data)
}

#' import mutation file
importMutationFileOnly <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0) {
  data <- read_csv(mutation_file, show_col_types = FALSE)
  data <- data %>% filter(alt_reads/total_reads>=vaf_thresh, alt_reads>=alt_reads_thresh)
  output_data <- list()
  
  output_data$y <- as.matrix(data[c("mutation", "sample", "alt_reads")] %>% pivot_wider(names_from = sample, values_from = alt_reads, values_fill = 0))
  rownames(output_data$y) <- output_data$y[,'mutation']
  output_data$y <- output_data$y[,-1, drop=FALSE]
  rowname = rownames(output_data$y)
  colname = colnames(output_data$y)
  output_data$y <- matrix(as.numeric(output_data$y), ncol = ncol(output_data$y))
  rownames(output_data$y) = rowname
  colnames(output_data$y) = colname
  
  output_data$MutID <- rowname
  
  output_data$n <- as.matrix(data[c("mutation", "sample", "total_reads")] %>% pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0))
  rownames(output_data$n) <- output_data$n[,'mutation']
  output_data$n <- output_data$n[,-1, drop=FALSE]
  rowname = rownames(output_data$n)
  colname = colnames(output_data$n)
  output_data$n <- matrix(as.numeric(output_data$n), ncol = ncol(output_data$n))
  rownames(output_data$n) = rowname
  colnames(output_data$n) = colname
  
  if (any((output_data$y - output_data$n) > 0)) {
    warning("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
    stop()
  }
  
  if (any(output_data$n==0)) {
    print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- round(mean(output_data$n))
  }
  
  output_data$icn <- as.matrix(data[c("mutation", "sample", "tumor_integer_copy_number")] %>% pivot_wider(names_from = sample, values_from = tumor_integer_copy_number, values_fill = 2))
  rownames(output_data$icn) <- output_data$icn[,'mutation']
  output_data$icn <- as.numeric(output_data$icn[,2])
  # rowname = rownames(output_data$icn)
  # colname = colnames(output_data$icn)
  # output_data$icn <- matrix(as.numeric(output_data$icn), ncol = ncol(output_data$icn))
  # rownames(output_data$icn) = rowname
  # colnames(output_data$icn) = colname
  

  
  # output_data$mtp <- estimateMultiplicityMatrix(output_data)[,1]
  
  # rowname = rownames(output_data$mtp)
  # colname = colnames(output_data$mtp)
  # output_data$mtp <- matrix(as.numeric(output_data$mtp), ncol = ncol(output_data$mtp))
  # rownames(output_data$mtp) = rowname
  # colnames(output_data$mtp) = colname
  
  output_data$cncf <- as.matrix(data[c("mutation", "sample", "cncf")] %>% pivot_wider(names_from = sample, values_from = cncf, values_fill = 1))
  rownames(output_data$cncf) <- output_data$cncf[,'mutation']
  output_data$cncf <- output_data$cncf[,-1,drop=FALSE]
  rowname = rownames(output_data$cncf)
  colname = colnames(output_data$cncf)
  output_data$cncf <- matrix(as.numeric(output_data$cncf), ncol = ncol(output_data$cncf))
  rownames(output_data$cncf) = rowname
  colnames(output_data$cncf) = colname
  
  output_data$S = ncol(output_data$y)
  output_data$I = nrow(output_data$y)
  
  output_data$tcn = output_data$icn * output_data$cncf + 2 * ( 1 - output_data$cncf)
  
  if ("major_integer_copy_number" %in% colnames(data)) {
    output_data$mtp <- as.matrix(data[c("mutation", "sample", "major_integer_copy_number")] %>% pivot_wider(names_from = sample, values_from = major_integer_copy_number, values_fill = 1))
    rownames(output_data$mtp) <- output_data$mtp[,'mutation']
    output_data$mtp <- as.numeric(output_data$mtp[,2])
  } else {
    output_data$mtp <- estimateMultiplicityMatrix(output_data)[,1]
  }
  
  if ("purity" %in% colnames(data)) {
    output_data$purity <- as.matrix(data[c("mutation", "sample", "purity")] %>% pivot_wider(names_from = sample, values_from = purity, values_fill = 0))
    rownames(output_data$purity) <- output_data$purity[,'mutation']
    output_data$purity <- output_data$purity[,-1, drop=FALSE]
    rowname = rownames(output_data$purity)
    colname = colnames(output_data$purity)
    output_data$purity <- matrix(as.numeric(output_data$purity), ncol = ncol(output_data$purity))
    rownames(output_data$purity) = rowname
    colnames(output_data$purity) = colname
    output_data$purity = colSums(output_data$purity) / colSums(!!output_data$purity)
  } else {
    output_data$purity = rep(0.8, ncol(output_data$y))
  }
  
  # mutation_position = unique(data[c("mutation", "chrom", "start", "end")])
  # if (!all(sort(unique(data[c("mutation", "chrom", "start", "end")]) %>% pull(mutation)) == sort(output_data$MutID))) {
  #   warning("Some mutations may have duplicated chromosome information; keeping the first occurence.")
  #   mutation_position = mutation_position[match(unique(mutation_position$mutation), mutation_position$mutation), ]
  # }
  # 
  # output_data$position = mutation_position
  
  
  return(output_data)
}

simulation_hidden <- function(mcf=0.9, icn=2, minor_cn=1, depth=30, num_SNV=30, seed=NULL, normal_prop=1) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  tcn = 2 * (1-mcf) + icn * mcf
  
  # if a proportion of SNVs on CNA is actually on CN-neutral region
  num_neutral = round(num_SNV * normal_prop)
  SNV_depth_neutral = rpois(n = num_neutral, lambda = depth)
  SNV_alt_neutral = numeric(num_neutral)
  for (i in seq_len(num_neutral)) {
    SNV_alt_neutral[i] <- rbinom(n=1, size=SNV_depth_neutral[i], prob=0.5)
  }
  
  vaf_neutral = SNV_alt_neutral / SNV_depth_neutral
  
  # Actual tumor proportion
  num_SNV = num_SNV - num_neutral
  
  # generate depth for each SNV
  depth_total = rpois(n = num_SNV, lambda = depth * tcn / 2)
  
  # assign each SNV to one copy
  SNV_assignment = sample(c(1, 2), size = num_SNV, replace = TRUE) # which segment
  
  # generate depth for germline
  SNV_depth_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_depth_germline[i] <- rbinom(n=1, size=depth_total[i], prob=2 * (1-mcf)/tcn)
  }
  
  SNV_alt_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_alt_germline[i] <- rbinom(n=1, size=SNV_depth_germline[i], prob=0.5)
  }
  
  vaf_germline <- SNV_alt_germline/SNV_depth_germline
  
  # generate depth for tumor
  SNV_depth_tumor = depth_total - SNV_depth_germline 
  SNV_alt_tumor = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    if (icn == 0) {
      SNV_alt_tumor[i] <- 0
    } else {
      tumor_vaf = minor_cn / icn
      if (SNV_assignment[i] == 1) {
        tumor_vaf = 1 - (minor_cn / icn)
      } 
      SNV_alt_tumor[i] <- rbinom(n=1, size=SNV_depth_tumor[i], prob=tumor_vaf)
    }
  }
  
  vaf_tumor = c((SNV_alt_germline + SNV_alt_tumor) / depth_total, vaf_neutral)
  
  is.unimodal(vaf_tumor)
}

