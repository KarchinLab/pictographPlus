#' sample presence for MCMC
#' @export
separateMutationsBySamplePresence <- function(input_data, min_mutation_per_cluster=5) {
  # returns list of lists -- 
  # each item of list contains input data for a mutation sample presence set 
  # original mutation indices from input_data are recorded in $mutation_indices
  pres <- ifelse(input_data$y > 0 & !input_data$is_cn, 1, 0) + ifelse(input_data$is_cn & input_data$tcn != 2, 1, 0)

  pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
  types <- sort(names(table(pat)), decreasing=TRUE)

  type_indices <- lapply(types, function(x) which(pat == x))
  names(type_indices) <- types
  
  sep_list <- list()
  toWarn = FALSE
  for (t in seq_len(length(types))) {
    sep_list[[types[t]]] <- list(pattern = types[t],
                                 mutation_indices = type_indices[[types[t]]],
                                 y = input_data$y[type_indices[[types[t]]], ,drop=FALSE],
                                 n = input_data$n[type_indices[[types[t]]], ,drop=FALSE],
                                 tcn = input_data$tcn[type_indices[[types[t]]], ,drop=FALSE],
                                 is_cn = input_data$is_cn[type_indices[[types[t]]]],
                                 cncf = input_data$cncf[type_indices[[types[t]]], ,drop=FALSE],
                                 mtp = input_data$mtp[type_indices[[types[t]]]],
                                 icn = input_data$icn[type_indices[[types[t]]]],
                                 MutID = input_data$MutID[type_indices[[types[t]]]],
                                 purity = input_data$purity)
    if (length(sep_list[[types[t]]]$MutID) < min_mutation_per_cluster) {toWarn = T}
  }
  
  if (toWarn) {warning("Using sample_presence mode, but exists clusters where number of mutations in the cluster is less than min_mutation_per_cluster. Please consider to set sample_presence = FALSE")}
  
  return(sep_list)
}