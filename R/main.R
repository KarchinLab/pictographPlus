#' main steps to run pictograph2
#' @export
runPictograph2 <- function(){
  # read in files
  input_data <- importFiles('./inst/extdata/sim_v2_snv.csv', './inst/extdata/sim_v2_cn.csv')

  # separate by sample presence
  sep_list <- separateBySamplePresence(input_data)

  # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  cna_est <- estimateCNA(input_data$tcn)

  # estimation of copy number cellular fraction
  cncf_est <- estimateCNCF(input_data$tcn, cna_est)
  
  # estimate multiplicity
  m_est <- estimateMultiplicity(input_data, cna_est, cncf_est)

  # estimation of mutation cellular fraction for each mutation
  mcf_est <- estimateMutation(input_data, sep_list, cna_est, cncf_est)

  # combined estimation
  total_est <- rbind(mcf_est, cncf_est)

  # separate mcf by sample presence
  total_sep_list <- separateMutationsByMCF(total_est)

  # return(input_data)
}
