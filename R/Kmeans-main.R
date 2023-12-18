#' main steps to run pictograph2
#' @export
kmeansMain <- function(){
  # read in files
  input_data <- importFiles('./inst/extdata/sim_v2_snv.csv', './inst/extdata/sim_v2_cn.csv')

  # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  cna_est <- estimateCNA(input_data)

  # estimation of copy number cellular fraction
  cncf_est <- estimateCNCF(input_data, cna_est)
  
  # estimate multiplicity
  # m_est <- estimateMultiplicity1(input_data, cna_est, cncf_est)

  # estimation of mutation cellular fraction for each mutation
  mcf_est <- estimateMutation(input_data, cna_est, cncf_est)

  # combine estimation for mcf and cncf
  total_est <- rbind(mcf_est, cncf_est)

  # separate mcf by sample presence
  total_sep_list <- separateMutationsByMCF(total_est)

  # test '101' box
  clusters = KmeansClustering(total_est, total_sep_list)
  
  clusterMCF = estimateClusterMCF(total_est, clusters)
  clusterAssignment = getClusterAssignment(total_est, clusters)
  
  # load_all("../pictograph")
  # w_mat <- clusterMCF
  # lineage_precedence_thresh=0.1
  # sum_filter_thresh=0.2
  # graph_G_pre <- prepareGraph(w_mat, lineage_precedence_thresh)
  # graph_G <- filterEdgesBasedOnCCFs(graph_G_pre, w_mat, thresh = lineage_precedence_thresh)
  # enumerateSpanningTreesModified(graph_G, w_mat, sum_filter_thresh = sum_filter_thresh)
  # scores <- calcTreeScores(chains$w_chain, all_spanning_trees)
  # best_tree <- all_spanning_trees[[which.max(scores)]]
  # plotTree(best_tree)
  
  # pattern = '101'
  # mutation_indices=sep_list$mutation[[pattern]]
  # tempBox <- list(pattern=pattern,
  #                 mutation_indices=mutation_indices,
  #                 y = input_data$y[mutation_indices, ,drop=FALSE],
  #                 n = input_data$n[mutation_indices, ,drop=FALSE],
  #                 cna = cna_est[mutation_indices, ,drop=FALSE],
  #                 cncf = cna_est[mutation_indices, ,drop=FALSE]
  #                 m = m_est[mutation_indices, ,drop=FALSE])

  # # return(input_data)
}


