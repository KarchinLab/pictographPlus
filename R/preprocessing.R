#' read input data file and store in PICTOGRAPH input format
#' @export
#' @param mutation_file mutation data file that contains columns "sample", "mutation", "chrom", "start", "end", total_reads", and "alt_reads";
#' @param copy_number_file copy number file that contains columns "sample", "chrom", "start", "end", "tcn"
importFiles <- function(mutation_file, copy_number_file=NULL, alt_reads_thresh = 0, vaf_thresh = 0, cnv_max_dist=2000, cnv_max_percent=0.10, tcn_normal_range=c(1.8, 2.2), smooth_cnv=F, autosome=T) {
  
  # keep mutations if alt_reads >= alt_reads_thresh and vaf >= vaf_thresh
  mutation_data = importMutationFile(mutation_file, alt_reads_thresh, vaf_thresh)
  
  copy_number_data = importCopyNumberFile(copy_number_file, cnv_max_dist, cnv_max_percent, tcn_normal_range, smooth_cnv, autosome)
  warning("importFiles: CNA not checked for overlap yet; user need to make sure CNA seperate")
  mutation_data$tcn = copy_number_data$tcn
  # mutation_data$minor = copy_number_data$minor
  # mutation_data$tcn = mutation_data$major + mutation_data$minor
  mutation_data$overlap = resolveOverlap(mutation_data)
  warning("resolveOverlap: need to check whether a mutation overlaps with two CNA segs")

  return(mutation_data)
}

#' check whether a mutation overlaps a CNA region
resolveOverlap <- function(mutation_data) {
  mut_count = nrow(mutation_data$position)
  cna_count = nrow(mutation_data$tcn)
  output = matrix(0, nrow=mut_count, ncol=cna_count)
  rownames(output) = mutation_data$position$mutation
  colnames(output) = rownames(mutation_data$tcn)
  # print(output)
  cnas = rownames(mutation_data$tcn)
  mutations = mutation_data$position
  for (i in seq_len(mut_count)) {
    for (j in seq_len(cna_count)) {
      chrSplit = strsplit(cnas[j], split='-')[[1]]
      chr = chrSplit[1]
      start_pos = strtoi(chrSplit[2])
      end_pos = strtoi(chrSplit[3])

      if (chr==mutations[i,2]) {
        if (strtoi(mutations[i,3]) >= start_pos && strtoi(mutations[i,4]) <= end_pos) {
          output[i,j] = 1
        }
      }
    }
  }
  return(output)
}

#' merge segments with similar coordinates
#' the maximum allowed distance between the start or end position of two segments is the max(cnv_max_dist, min(length of two segments)*cnv_max_percent)
smoothCNV <- function(data, cnv_max_dist=2000, cnv_max_percent=0.10) {
  # smooth_data <- as_tibble(data.frame(matrix(nrow=0, ncol=ncol(data))))
  # colnames(smooth_data) <- colnames(data)
  
  # samples <- unique(data$sample)
  data <- data %>% arrange(chrom, start)
  indexList = seq_len(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    if (i %in% indexList) {
      # print(data[i,])
      max_dis = max(cnv_max_dist, cnv_max_percent*(data[i,]$end-data[i,]$start))
      index_selected = which((data$chrom==data[i,]$chrom)&(abs(data$start-data[i,]$start)<max_dis)&(abs(data$end-data[i,]$end)<max_dis))
      
      # list of cnv segments able to merge
      cnv_selected = data[index_selected,]
      # start = min(cnv_selected$start)
      # end = max(cnv_selected$end)
      data[index_selected,]$start = min(cnv_selected$start)
      data[index_selected,]$end = max(cnv_selected$end)
      indexList <- indexList[!(indexList %in% index_selected)]
    }
  }
  warning("Add split CNV function in smoothCNV under preprocessing.R")
  return(data)
}

#' import copy number file
#' @param cnv_max_dist: maximum of distance allowed between two segments to assign as the same one
#' @param cnv_max_percent: maximum percentage of distance allowed between two segments to assign as the same one
#' @param smooth_cnv: process input CNV to merge  segments with similar distance
importCopyNumberFile <- function(copy_number_file, cnv_max_dist=2000, cnv_max_percent=0.10, tcn_normal_range=c(1.8, 2.2), smooth_cnv=T, autosome=T) {
  data <- read_csv(copy_number_file, show_col_types = FALSE)
  # data <- data %>% arrange(chrom, start)
  # SMOOTH SEGMENTS
  if (smooth_cnv) {
    data <- smoothCNV(data, cnv_max_dist, cnv_max_percent)
  }
  
  data$CNA = paste(data$chrom, data$start, data$end, sep = '-')
  
  # keep only autosome
  if (autosome) {
    data <- data[!(data$chrom %in% c('chrX', 'chrY', 'chrM')),]
  }
  
  output_data <- list()
  output_data$tcn = as.matrix(data[c("sample", "CNA", "tcn")] %>% pivot_wider(names_from = sample, values_from = tcn, values_fill = 2))
  
  
  warning("FILL MISSING TCN WITH 2; NEED TO FIX FOR CHRX/Y FOR VALUES_FILL in importCopyNumberFile in preprocessing.R")
  rownames(output_data$tcn) <- output_data$tcn[,'CNA']
  output_data$tcn <- output_data$tcn[,-1, drop=FALSE]
  rowname = rownames(output_data$tcn)
  colname = colnames(output_data$tcn)
  output_data$tcn <- matrix(as.numeric(output_data$tcn), ncol = ncol(output_data$tcn))
  to_keep_index = which(rowSums(output_data$tcn<tcn_normal_range[1]|output_data$tcn>tcn_normal_range[2])>0)
  output_data$tcn = output_data$tcn[to_keep_index,]
  rownames(output_data$tcn) = rowname[to_keep_index]
  colnames(output_data$tcn) = colname

  # output_data$minor = as.matrix(data[c("sample", "CNA", "minor")] %>% pivot_wider(names_from = sample, values_from = minor, values_fill = 1))
  # rownames(output_data$minor) <- output_data$minor[,'CNA']
  # output_data$minor <- output_data$minor[,-1, drop=FALSE]
  # rowname = rownames(output_data$minor)
  # colname = colnames(output_data$minor)
  # output_data$minor <- matrix(as.numeric(output_data$minor), ncol = ncol(output_data$minor))
  # rownames(output_data$minor) = rowname
  # colnames(output_data$minor) = colname

  return(output_data)
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

  # output_data$y[output_data$y / output_data$n < vaf_thresh] = 0
  # output_data$y[output_data$y < alt_reads_thresh] = 0

  mutation_position = unique(data[c("mutation", "chrom", "start", "end")])
  if (!all(sort(unique(data[c("mutation", "chrom", "start", "end")]) %>% pull(mutation)) == sort(output_data$MutID))) {
    warning("Some mutations may have duplicated chromosome information; keeping the first occurence.")
    mutation_position = mutation_position[match(unique(mutation_position$mutation), mutation_position$mutation), ]
  }

  output_data$position = mutation_position
  return(output_data)
}

