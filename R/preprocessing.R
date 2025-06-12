#' data processing for mutation and copy number files
#' 
#' read input data file and store in required format
#' 
#' @export
#' @param mutation_file a csv file that include information for SSMs. See vignette for details.
#' @param copy_number_file a csv file that include information for CNA. See vignette for details.
#' @param SNV_file a csv file that include information for germline heterozygous SNVs. See vignette for details.
#' @param outputDir output directory for saving all files.
#' @param LOH whether or not to include copy number segments that are copy neutral but LOH
#' @param purity_min minimum purity for tumor samples
#' @param alt_reads_thresh minimum number of alternative read count for a SSM to be included in the analysis
#' @param vaf_thresh minimum VAF for a SSM to be included in the analysis
#' @param cnv_min_length minimum length of copy number alterations for it to be included in analysis
#' @param tcn_normal_range range of total copy number considered as copy-neutral
#' @param filter_cnv whether or not to filter copy number alterations
#' @param smooth_cnv whether or not to process copy number alterations across samples to unify the segment start and end positions
#' @param autosome to only include autosomes
#' @param pval placeholder

importFiles <- function(mutation_file, 
                        copy_number_file=NULL, 
                        SNV_file=NULL, 
                        outputDir=NULL,
                        LOH = FALSE,
                        purity_min = 0.2,
                        alt_reads_thresh = 0,
                        vaf_thresh = 0,
                        cnv_min_length=1000000,
                        tcn_normal_range=c(1.75, 2.3),
                        filter_cnv = T,
                        smooth_cnv= T,
                        autosome=T,
                        pval=0.05,
                        depth=NULL
                        ) {
  
  # set output directory to current directory if outputDir is NULL
  if (is.null(outputDir)) {
    outputDir = getwd()
  }
  
  if (!is.null(copy_number_file)) {
    # keep mutations if alt_reads >= alt_reads_thresh and vaf >= vaf_thresh
    mutation_data = importMutationFile(mutation_file, alt_reads_thresh, vaf_thresh, purity_min, depth)
    
    # sample orders
    name_order <- colnames(mutation_data$y)
    
    copy_number_data = importCopyNumberFile(copy_number_file, 
                                            outputDir, 
                                            SNV_file, 
                                            LOH, 
                                            name_order,
                                            cnv_min_length, 
                                            tcn_normal_range,
                                            filter_cnv,
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
    mutation_data = importMutationFileOnly(mutation_file, alt_reads_thresh, vaf_thresh, purity_min, depth)
    mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)))
    mutation_data$cnnull <- TRUE
  }
  
  mutation_data$y[mutation_data$is_cn==1,] <- pmax(mutation_data$y[mutation_data$is_cn==1,], mutation_data$n[mutation_data$is_cn==1,]-mutation_data$y[mutation_data$is_cn==1,])
  
  return(mutation_data)
}

#' check whether a mutation overlaps a CNA region
resolveOverlap <- function(mutation_data) {
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

resegmentation <- function(data, cnv_min_length=2000) {
  breakpoints <- data %>%
    select(chrom, start, end) %>%
    pivot_longer(cols = c(start, end), names_to = "position_type", values_to = "position") %>%
    distinct(chrom, position) %>%
    arrange(chrom, position)
  
  segments <- breakpoints %>%
    group_by(chrom) %>%
    arrange(position) %>%
    mutate(next_position = lead(position)) %>%
    filter(!is.na(next_position)) %>%
    select(chrom, start.x = position, end.x = next_position)
  
  expanded_data <- segments %>%
    left_join(data, by = "chrom") %>%
    filter(start.x >= start & end.x <= end) %>%
    select(-start, -end) %>%
    rename(start=start.x, end=end.x) %>%
    filter(end - start > cnv_min_length)%>%
    arrange(sample, chrom, start)
  
  return(expanded_data)
}

#' @import GenomicRanges IRanges dplyr tidyr
reSegCNV <- function(data, cnv_min_length) {
  
  # data_merged <- data %>%
  #   # 1) sort
  #   arrange(sample, chrom, start) %>%
  #   
  #   # 2) flag mergeable to previous
  #   group_by(sample, chrom) %>%
  #   mutate(
  #     tcn_diff = abs(tcn - lag(tcn)),
  #     baf_diff = abs(baf - lag(baf)),
  #     mergeable = (tcn_diff <= 0.2 & baf_diff <= 0.1) %>% replace_na(FALSE),
  #     
  #     # 3) assign a new group whenever NOT mergeable
  #     grp = cumsum(!mergeable),
  #     
  #     # compute segment length for weighting
  #     seg_length = end - start
  #   ) %>%
  #   
  #   # 4) collapse each run‐group
  #   group_by(sample, chrom, grp) %>%
  #   summarise(
  #     start = min(start),
  #     end   = max(end),
  #     tcn   = weighted.mean(tcn, w = seg_length),
  #     baf   = weighted.mean(baf, w = seg_length),
  #     .groups = "drop"
  #   ) %>%
  #   arrange(sample, chrom, start)
  
  grl <- split(data, data$sample) |> 
    lapply(function(df) {
      GRanges(
        seqnames = df$chrom,
        ranges   = IRanges(start = df$start, end = df$end),
        tcn      = df$tcn,
        baf      = df$baf
      )
    })
  grl <- GRangesList(grl)
  
  # — 2) Create the consensus “tiny bins”
  all_gr        <- unlist(grl, use.names = FALSE)
  consensus_gr  <- disjoin(all_gr)
  
  # — 3) Annotate each bin with each sample’s tcn & baf
  for (s in names(grl)) {
    hits <- findOverlaps(consensus_gr, grl[[s]])
    mcols(consensus_gr)[[paste0("tcn_", s)]] <- NA_real_
    mcols(consensus_gr)[[paste0("baf_", s)]] <- NA_real_
    mcols(consensus_gr)[[paste0("tcn_", s)]][queryHits(hits)] <- mcols(grl[[s]])$tcn[subjectHits(hits)]
    mcols(consensus_gr)[[paste0("baf_", s)]][queryHits(hits)] <- mcols(grl[[s]])$baf[subjectHits(hits)]
  }
  
  # — 4) Turn into a tidy data.frame
  df <- data.frame(
    seqnames = as.character(seqnames(consensus_gr)),
    start    = start(consensus_gr),
    end      = end(consensus_gr),
    mcols(consensus_gr),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # — 5) Prepare for breakpoint logic
  tcn_cols <- grep("^tcn_", names(df), value = TRUE)
  tol      <- 0.4   # adjust: 0 for exact, ~0.3–0.5 to smooth small jumps
  
  df2 <- df %>%
    mutate(
      across(all_of(tcn_cols),   ~ coalesce(., 2)),
      across(starts_with("baf_"), ~ coalesce(., 0.5))
    ) %>%
    mutate(width = end - start + 1) %>%
    group_by(seqnames) %>%
    mutate(
      # pack TCN into a matrix for diff’ing
      tcn_mat    = list(as.matrix(across(all_of(tcn_cols)))),
      break_flag = {
        m     <- tcn_mat[[1]]
        diffs <- abs(m[-1L, , drop=FALSE] - m[-nrow(m), , drop=FALSE])
        c(TRUE, apply(diffs, 1, function(x) any(x > tol)))
      },
      seg_id = cumsum(break_flag)
    ) %>%
    ungroup() %>%
    select(-tcn_mat, -break_flag)
  
  # — 6) Summarise into shared segments and re-pivot, using names_pattern
  final <- df2 %>%
    group_by(seqnames, seg_id) %>%
    summarise(
      start = min(start),
      end   = max(end),
      across(all_of(tcn_cols),  ~ sum(. * width, na.rm=TRUE) / sum(width, na.rm=TRUE)),
      across(starts_with("baf_"), ~ sum(. * width, na.rm=TRUE) / sum(width, na.rm=TRUE)),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols       = matches("^(tcn|baf)_"),
      names_to   = c("metric","sample"),
      names_pattern = "^(tcn|baf)_(.+)$",
      values_to  = "value"
    ) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    rename(chrom = seqnames) %>%
    select(sample, chrom, start, end, tcn, baf) %>%
    filter(end - start > cnv_min_length)
  
  return(final)
}

#' import copy number file
#' @import GenomicRanges IRanges dplyr tidyr
importCopyNumberFile <- function(copy_number_file, outputDir, SNV_file=NULL, 
                                 LOH=FALSE, name_order=NULL, cnv_min_length=1000000, 
                                 tcn_normal_range=c(1.75, 2.3), filter_cnv = T, 
                                 smooth_cnv=F, autosome=T, pval=0.05) {
  
  data <- read_csv(copy_number_file, show_col_types = FALSE) # read copy number csv file
  
  if (nrow(data)==0) {
    return(NULL)
  }
  
  if (smooth_cnv) {
    if ("baf" %in% colnames(data)) {
      data <- reSegCNV(data, cnv_min_length)
    } else {
      data <- resegmentation(data, cnv_min_length)
    }
  }
  
  if ("baf" %in% colnames(data)) {
    message("inferring allele-specific copy number using BAF")
    
    if (filter_cnv) {
      
      data <- data %>%
        group_by(chrom, start, end) %>%
        filter(
          sum(
            tcn < tcn_normal_range[1] |
              tcn > tcn_normal_range[2] |
              baf < 0.3 |
              baf > 0.7
          ) >= 2
        ) %>%
        ungroup()
      
      data <- data %>% 
        filter(!(tcn >= 1.9 & tcn <= 2.1 & baf >= 0.4 & baf <= 0.6))
    }
    
  } else if (!is.null(SNV_file)) {
    message("inferring allele-specific copy number using heterozygous SNVs")
    
    # check unimodality in both normal and tumor sample
    data <- check_sample_CNA(data, outputDir, SNV_file, LOH, tcn_normal_range=tcn_normal_range, pval=pval) 
    if (filter_cnv) {
      data <- data[data$to_keep==1,] # keep rows is to_keep is 1
    }
  } else {
    stop("Please provide either a baf column in the copy number file or supply a SNV file with heterozygous SNV counts")
  }
  
  data$tcn[data$tcn==2] <- 2.001 # make tcn=2 -> 2.01 to avoid confusion during sample presence
  
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
    
    tcn_tot <- matrix(500, nrow(output_data), ncol(output_data))
    rownames(tcn_tot) <- rownames(output_data)
    colnames(tcn_tot) <- colnames(output_data)
    
    tcn_tot <- add_missing_column(name_order, tcn_tot, 500)
    
    tcn_alt <- matrix(round(tcn_tot * baf), nrow(output_data),ncol(output_data))
    tcn_alt <- pmax(tcn_alt, tcn_tot - tcn_alt)
    rownames(tcn_alt) <- rownames(output_data)
    colnames(tcn_alt) <- colnames(output_data)
    
    # tcn_alt <- add_missing_column(name_order, tcn_alt, 500)
    
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
    
    tcn_tot <- round(tcn_tot) * 5
    tcn_alt <- round(tcn_alt) * 5
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

#' copy number quality check
#' 
#' Check if a copy-neutral segment is a LOH event by checking the distribution of germline heterozygous mutations
#' @import LaplacesDemon parallel diptest
#' 
check_sample_CNA <- function(data, outputDir, SNV_file, LOH, tcn_normal_range=c(1.75, 2.3), pval=0.05) {

  SNV_data <- read_csv(SNV_file)
  
  if ("chroms" %in% colnames(SNV_data)) {
    SNV_data <- SNV_data %>%
      rename(chrom = chroms)
  }
  
  samples <- unique(data$sample)

  pdf(paste(outputDir, "sample_modality.pdf", sep="/"), width = 12, height = 18)
  par(mfrow=c(6,3), cex.main = 1, cex.axis = 1, cex.lab = 1, cex = 1)
  
  rowLen = nrow(data)

  tcn_ref <- rep(0, nrow(data))
  tcn_alt <- rep(0, nrow(data))

  to_keep_index = c()

  for (i in 1:rowLen) {
    toSkip = FALSE
    
    SNV_temp <- SNV_data %>% filter(chrom==data[i,]$chrom & position>=data[i,]$start & position<=data[i,]$end)
    sample <- data[i,]$sample
    
    if (nrow(SNV_temp) > 20) {
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
            if (tumor_test$p.value == 0) {
              if (LOH) {
                to_keep_index <- c(to_keep_index, 1)
              } else {
                to_keep_index <- c(to_keep_index, 0)
              }
              plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, ", pval: ", tumor_test$p.value, sep=""))

            } else { # disregard if only one peak
              # plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, ", pval: ", tumor_test$p.value, sep=""))
              to_keep_index <- c(to_keep_index, 0)
            }
          } else { # keep CNA because tcn not in normal range
            # if (tumor_test$p.value < pval/nrow(data)) {
            pval1 = 0.0001
            if (is.bimodal(vaf)) {
              pval1 = 0.001
            }
            if (tumor_test$p.value < pval1) {
              to_keep_index <- c(to_keep_index, 1)
              plot(density(vaf), xlim=c(0,1), main = paste(data[i,]$sample, "\n", data[i,]$chrom, ":", data[i,]$start, "-", data[i,]$end, "\n tcn: ", data[i,]$tcn, ", pval: ", tumor_test$p.value, sep=""))
            } else {
              if (data[i,]$tcn > 3 || data[i,]$tcn < 1.4) {
                to_keep_index <- c(to_keep_index, 1)
              } else {
                to_keep_index <- c(to_keep_index, 0)
              }
            }
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
importMutationFile <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0, purity_min=0.2, depth=300) {
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
    output_data$purity <- pmax(output_data$purity, purity_min)
  } else {
    output_data$purity = rep(0.8, ncol(output_data$y))
  }
  
  if (any((output_data$y - output_data$n) > 0)) {
    stop("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
  }

  if (any(output_data$n==0)) {
    # print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    # output_data$n[output_data$n==0] <- round(mean(output_data$n[output_data$n!=0]))
    if (is.null(depth)) {
      output_data$n <- update_n(output_data$n)
    } else {
      output_data$n[output_data$n==0] <- depth
    }
    
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

#' import only mutation file
importMutationFileOnly <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0, purity_min=0.2, depth=300) {
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
    # print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    # output_data$n[output_data$n==0] <- 300
    # output_data$n <- update_n(output_data$n)
    if (is.null(depth)) {
      output_data$n <- update_n(output_data$n)
    } else {
      output_data$n[output_data$n==0] <- depth
    }
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
    output_data$purity <- pmax(output_data$purity, purity_min)
  } else {
    output_data$purity = rep(0.8, ncol(output_data$y))
  }
  
  return(output_data)
}

update_n <- function(m) {
  
  int_means <- round(rowSums(m) / rowSums(m != 0))
  
  # 2. Replace zeros in each row with the corresponding int_mean
  m2 <- m  # copy
  for(i in seq_len(nrow(m2))) {
    zero_locs <- m2[i, ] == 0
    m2[i, zero_locs] <- int_means[i]
  }
  
  return(m2)
}
