#' read input data file and store in required format
#' 
#' @export
#' @param mutation_file mutation data file; see inst/extdata/examples/*_snv.csv for examples
#' @param copy_number_file copy number file; see inst/extdata/examples/*_cna.csv for examples
#' @param SNV_file SNV file for germline heterozygous SNVs; see inst/extdata/examples/*_SNP.csv for examples
#' @param outputDir output directory for saving output data
importFiles <- function(mutation_file, 
                        copy_number_file=NULL, 
                        SNV_file=NULL, 
                        outputDir=NULL,
                        alt_reads_thresh = 0, # to be tested
                        vaf_thresh = 0, # to be tested
                        cnv_max_dist=2000, # to be tested
                        cnv_max_percent=0.30, # to be tested
                        tcn_normal_range=c(1.8, 2.2), # to be tested
                        smooth_cnv=F, # to be tested
                        autosome=T, # to be tested
                        pval=0.05 # to be tested
                        ) {
  
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
                                            name_order,
                                            cnv_max_dist, 
                                            cnv_max_percent, 
                                            tcn_normal_range, 
                                            smooth_cnv, 
                                            autosome, 
                                            pval)
    
    if (is.null(copy_number_data)) {
      
      mutation_data$icn <- rep(2, nrow(mutation_data$y))

      mutation_data$mtp <- rep(1, nrow(mutation_data$y))

      mutation_data$cncf <- mutation_data$y
      mutation_data$cncf[] <- 0
      
      mutation_data$tcn = mutation_data$icn * mutation_data$cncf + 2 * ( 1 - mutation_data$cncf)
      
      mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)))
      
      mutation_data$cnnull <- TRUE
      
    } else {
      mutation_data$tcn = copy_number_data$tcn[, name_order, drop=FALSE]
      
      # bind SSM and CNA
      mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)), rep(1,nrow(copy_number_data$tcn)))
      mutation_data$y <- rbind(mutation_data$y, copy_number_data$tcn_alt[, name_order, drop=FALSE])
      mutation_data$n <- rbind(mutation_data$n, copy_number_data$tcn_tot[, name_order, drop=FALSE])
      mutation_data$MutID <- c(mutation_data$MutID, rownames(copy_number_data$tcn))
      mutation_data$I <- mutation_data$I + nrow(copy_number_data$tcn)
      
      mutation_data$overlap = resolveOverlap(mutation_data)
      mutation_data$tcn <- mutation_data$overlap %*% mutation_data$tcn
      # mutation_data$tcn[mutation_data$tcn==0] <- 2.0
      
      overlap <- mutation_data$overlap
      colnames(overlap) <- sapply(colnames(mutation_data$overlap), function(col) {which(rownames(mutation_data$overlap)==col)})
      
      q <- vector("numeric", nrow(overlap))
      for (i in 1:nrow(overlap)) {
        q[i] <- ifelse(length(which(overlap[i,] == 1)) > 0, as.numeric(names(which(overlap[i,] == 1))[1]),i)
      }
      mutation_data$q <- q
      
      mutation_data$cnnull <- FALSE
    }
    
  } else {
    mutation_data = importMutationFileOnly(mutation_file, alt_reads_thresh, vaf_thresh)
    mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)))
    mutation_data$cnnull <- TRUE
  }
  return(mutation_data)
}

#' check whether a mutation overlaps a CNA region
resolveOverlap <- function(mutation_data) {
  # NOTE: TO DO
  # warning("resolveOverlap: need to check whether a mutation overlaps with two CNA segs")
  mut_count = nrow(mutation_data$y)
  cna_count = nrow(mutation_data$tcn)
  output = matrix(0, nrow=mut_count, ncol=cna_count)
  rownames(output) = rownames(mutation_data$y)
  colnames(output) = rownames(mutation_data$tcn)
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
#' 
#' the maximum allowed distance between the start or end position of two segments 
#' is the max(cnv_max_dist, min(length of two segments)*cnv_max_percent)
smoothCNV <- function(data, cnv_max_dist=2000, cnv_max_percent=0.30) {
  # NOTE: TO DO
  # message("improvements needed for smoothCNV under preprocessing.R; 
  # currently detecting two segments as the same if overlap and star/end position 
  # within certain distance; merge the two segments as the outcome")
  data <- data %>% 
    mutate(numeric_chrom = as.numeric(sub("chr", "", chrom))) %>% 
    arrange(numeric_chrom, start) %>%
    select(-numeric_chrom)
  
  indexList = seq_len(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    # print(i)
    if (i %in% indexList) {
      max_dis = max(cnv_max_dist, cnv_max_percent*(data[i,]$end-data[i,]$start))
      index_selected = which((data$chrom==data[i,]$chrom)&(data$start<=data[i,]$end)&(data[i,]$start<=data$end)&(abs(data$start-data[i,]$start)<max_dis)&(abs(data$end-data[i,]$end)<max_dis))
      if (length(index_selected) > 1) {
        # list of cnv segments able to merge
        cnv_selected = data[index_selected,]
        data[index_selected,]$start = min(cnv_selected$start)
        data[index_selected,]$end = max(cnv_selected$end)
        # data[index_selected,]$drivers = paste(unique(unlist(strsplit(cnv_selected$drivers, ";"))), collapse=";")
        # data[index_selected,]$genes = paste(unique(unlist(strsplit(cnv_selected$genes, ";"))), collapse=";")
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
importCopyNumberFile <- function(copy_number_file, outputDir, SNV_file=NULL, name_order=NULL, cnv_max_dist=2000, cnv_max_percent=0.30, tcn_normal_range=c(1.8, 2.2), smooth_cnv=F, autosome=T, pval=0.05) {
  
  data <- read_csv(copy_number_file, show_col_types = FALSE) # read copy number csv file
  
  if (nrow(data)==0) {
    return(NULL)
  }
  
  if ("baf" %in% colnames(data)) {
    message("inferring allele-specific copy number using BAF")
    
    data <- data %>% 
      filter(!(tcn >= tcn_normal_range[1] & tcn <= tcn_normal_range[2] & baf >= 0.3 & baf <= 0.7))
    
  } else if (!is.null(SNV_file)) {
    message("inferring allele-specific copy number using heterozygous SNVs")
    # check unimodality at both normal and tumor sample
    
    # NOTE: TO DO; Not working yet thus only incuding cna in certain range, but not LOH
    data <- check_sample_LOH(data, outputDir, SNV_file, tcn_normal_range=tcn_normal_range, pval=pval) 
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
    
    baf <- add_missing_column(name_order, baf, 0.5)
    
    tcn_tot <- matrix(2000, nrow(output_data), ncol(output_data))
    rownames(tcn_tot) <- rownames(output_data)
    colnames(tcn_tot) <- colnames(output_data)
    
    tcn_tot <- add_missing_column(name_order, tcn_tot, 2000)
    
    tcn_alt <- matrix(round(tcn_tot * baf), nrow(output_data),ncol(output_data))
    tcn_alt <- pmax(tcn_alt, tcn_tot - tcn_alt)
    rownames(tcn_alt) <- rownames(output_data)
    colnames(tcn_alt) <- colnames(output_data)
    
    # tcn_alt <- add_missing_column(name_order, tcn_alt, 1000)
    
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
    alt_mean <- round(mean(tcn_tot[tcn_tot!=0])/2)
    tot_mean <- round(mean(tcn_tot[tcn_tot!=0]))
    tcn_alt[tcn_tot==0] <- alt_mean
    tcn_tot[tcn_tot==0] <- tot_mean
    
    tcn_tot <- add_missing_column(name_order, tcn_tot, tot_mean)
    tcn_alt <- add_missing_column(name_order, tcn_alt, alt_mean)
    
    tcn_tot <- round(tcn_tot) * 10
    tcn_alt <- round(tcn_alt) * 10
  }
  
  return_data <- list()
  return_data$tcn <- output_data
  return_data$tcn_tot <- tcn_tot
  return_data$tcn_alt <- tcn_alt
  return(return_data)
}

#' add back samples with no CNA events 
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

#' Check if a copy-neutral segment is a LOH event by checking the distribution of germline heterozygous mutations
#' @import LaplacesDemon
#' @import parallel
#' @import diptest
#' 
check_sample_LOH <- function(data, outputDir, SNV_file, tcn_normal_range=c(1.5, 2.5), pval=0.05) {

  SNV_data <- read_csv(SNV_file)
  samples <- unique(data$sample)

  pdf(paste(outputDir, "sample_modality.pdf", sep="/"), width = 12, height = 18)
  par(mfrow=c(6,3), cex.main = 1, cex.axis = 1, cex.lab = 1, cex = 1)
  
  rowLen = nrow(data)

  tcn_ref <- rep(0, nrow(data))
  tcn_alt <- rep(0, nrow(data))

  to_keep_index = c()

  for (i in 1:rowLen) {
    toSkip = FALSE
    
    SNV_temp <- SNV_data %>% filter(chroms==data[i,]$chrom & position>=data[i,]$start & position<=data[i,]$end)
    sample <- data[i,]$sample
    
    if (nrow(SNV_temp) > 2) {
      vaf_germline <- SNV_temp$germline_alt / (SNV_temp$germline_alt + SNV_temp$germline_ref)
      germline_test <- dip.test(vaf_germline)
      alt = paste(sample, "alt", sep="_")
      ref = paste(sample, "ref", sep="_")
      
      vaf <- SNV_temp[[alt]] / (SNV_temp[[alt]] + SNV_temp[[ref]])
      
      tumor_test <- dip.test(vaf)
      
      if (is.unimodal(vaf)) { # if similar to germline distribution, alt and ref will be the sum of all
        tcn_alt[i] = mean(SNV_temp[[alt]])
        tcn_ref[i] = mean(SNV_temp[[ref]])
      } else if (is.trimodal(vaf)) { # else, cluster the vaf into two and take one cluster
        kmeans_result <- kmeans(vaf, centers = 3)
        cluster_number = which.max(kmeans_result$centers)
        tcn_alt[i] = mean(SNV_temp[[alt]][kmeans_result$cluster == cluster_number])
        tcn_ref[i] = mean(SNV_temp[[ref]][kmeans_result$cluster == cluster_number])
      } else if (is.bimodal(vaf)) {
        kmeans_result <- kmeans(vaf, centers = 2)
        cluster_number = which.max(kmeans_result$centers)
        tcn_alt[i] = mean(SNV_temp[[alt]][kmeans_result$cluster == cluster_number])
        tcn_ref[i] = mean(SNV_temp[[ref]][kmeans_result$cluster == cluster_number])
      } else {
        tcn_alt[i] = -1
        tcn_ref[i] = -1
        toSkip = TRUE
      }

      # check if a segment is within tcn_normal_range for LOH
      if (toSkip) {
        to_keep_index <- c(to_keep_index, 0)
      } else {
        if (!germline_test$p.value < pval/nrow(data)) { # proceed is germline distribution is unimodality
          if (data[i,]$tcn > tcn_normal_range[1] & data[i,]$tcn < tcn_normal_range[2]) { # check tcn within normal range
            if (tumor_test$p.value < pval/nrow(data)) {
              # NOTE: TO DO
              warning("check_sample_LOH diptest seems not working; not detecting LOH; more testing needed")
              # to_keep_index <- c(to_keep_index, 1)
              to_keep_index <- c(to_keep_index, 0)
              plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, ", pval: ", tumor_test$p.value, sep=""))

            } else { # disregard if only one peak
              to_keep_index <- c(to_keep_index, 0)
            }
          } else { # keep CNA because tcn not in normal range
            to_keep_index <- c(to_keep_index, 1)
            plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, sep=""))
          }
        } else { # proceed is germline distribution is not unimodality
          to_keep_index <- c(to_keep_index, 0)
        }
      }
    } else { # not enough germline position; disregard
      to_keep_index <- c(to_keep_index, 0)
    }
  }
  dev.off()
  
  data <- data %>% add_column(to_keep = to_keep_index, tcn_ref = tcn_ref, tcn_alt = tcn_alt)
  
  return(data)
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
    stop("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
  }

  if (any(output_data$n==0)) {
    # print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    # output_data$n[output_data$n==0] <- round(mean(output_data$n[output_data$n!=0]))
    output_data$n[output_data$n==0] <- 2000
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
    stop("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
  }
  
  if (any(output_data$n==0)) {
    print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- 2000
  }
  
  output_data$icn <- as.matrix(data[c("mutation", "sample", "tumor_integer_copy_number")] %>% pivot_wider(names_from = sample, values_from = tumor_integer_copy_number, values_fill = 2))
  rownames(output_data$icn) <- output_data$icn[,'mutation']
  output_data$icn <- output_data$icn[,2:ncol(output_data$icn)]
  output_data$icn <- matrix(as.numeric(output_data$icn), ncol=ncol(output_data$y))
  output_data$icn <- apply(output_data$icn, 1, function(x) {
    if (all(x == 2)) {
      return(2)
    } else {
      return(mean(x[x != 2]))
    }
  })

  output_data$cncf <- as.matrix(data[c("mutation", "sample", "cncf")] %>% pivot_wider(names_from = sample, values_from = cncf, values_fill = 0))
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
    output_data$mtp <- output_data$mtp[,2:ncol(output_data$mtp)]
    output_data$mtp <- matrix(as.numeric(output_data$mtp), ncol=ncol(output_data$y))
    output_data$mtp <- apply(output_data$mtp, 1, function(x) {
      if (all(x == 1)) {
        return(1)
      } else {
        return(mean(x[x != 1]))
      }
    })
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
  
  return(output_data)
}
