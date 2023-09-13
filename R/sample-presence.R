
#' separate mutations by sample presence using mcf
#' @export
separateMutationsByMCF <- function(mcf) {
  # returns list of lists --
  # each item of list contains mcf for a mutation sample presence set
  # original mutation indices from input_data are recorded in $mutation_indices
  pres <- ifelse(mcf > 0, 1, 0)
  pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
  types <- sort(names(table(pat)), decreasing=TRUE)
  if (length(types) == 1) {
    type_indices <- list()
    type_indices[[types]] <- seq_len(length(nrow(mcf_est)))
  } else {
    type_indices <- lapply(types, function(x) which(pat == x))
    names(type_indices) <- types
  }
  warning('consider merge one-sample cluster with other clusters')
  return(type_indices)
}

#' sample presence for MCMC
#' @export
separateMutationsBySamplePresence <- function(input_data, min_mutation_per_box=2) {
  # returns list of lists -- 
  # each item of list contains input data for a mutation sample presence set 
  # original mutation indices from input_data are recorded in $mutation_indices
  pres <- ifelse(input_data$y > 0, 1, 0)
  pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
  types <- sort(names(table(pat)), decreasing=TRUE)
  if (length(types) == 1) {
    type_indices <- list()
    type_indices[[types]] <- seq_len(nrow(input_data$y))
  } else {
    type_indices <- lapply(types, function(x) which(pat == x))
    names(type_indices) <- types
  }
  
  sep_list <- list()
  for (t in seq_len(length(types))) {
    sep_list[[types[t]]] <- list(pattern = types[t],
                                 mutation_indices = type_indices[[types[t]]],
                                 cncf = input_data$cncf[type_indices[[types[t]]], ,drop=FALSE],
                                 y = input_data$y[type_indices[[types[t]]], ,drop=FALSE],
                                 n = input_data$n[type_indices[[types[t]]], ,drop=FALSE],
                                 tcn = input_data$tcn[type_indices[[types[t]]], ,drop=FALSE],
                                 m = input_data$m[type_indices[[types[t]]], ,drop=FALSE])
    if (ncol(input_data$y) == 1) {
      break
    }
  }
  return(sep_list)
}


#' #' separate mutation and copy number data based on sample presence profile
#' #' export
#' separateBySamplePresence <- function(input_data) {
#' 
#'   mutationTypes = separateMutationsBySamplePresence(input_data)
#'   # print(mutationTypes)
#' 
#'   cnaTypes = separateCNAsBySamplePresence(input_data)
#'   # print(cnaTypes)
#' 
#'   totalTypes = combineMutationCNASamplePresence(input_data)
#' 
#'   sep_list <- list()
#'   sep_list$mutation <- mutationTypes
#'   sep_list$cna <- cnaTypes
#'   sep_list$total <- totalTypes
#'   print(sep_list)
#'   warning("need to test if number of type = 1")
#'   return(sep_list)
#' }


#' combine mutation and cna sample presence
# combineMutationCNASamplePresence <- function(input_data) {
#   mut_pres <- ifelse(input_data$y > 0, 1, 0)
#   cna_pres <- ifelse(input_data$tcn == 2, 0, 1)
#   pres <- rbind(mut_pres, cna_pres)
#   pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
#   types <- sort(names(table(pat)), decreasing=TRUE)
#   if (length(types) == 1) {
#     type_indices <- list()
#     type_indices[[types]] <- seq_len(nrow(pres))
#   } else {
#     type_indices <- lapply(types, function(x) which(pat == x))
#     names(type_indices) <- types
#   }
#   return(type_indices)
# }

#' separate mutations by sample presence
# separateMutationsBySamplePresence <- function(input_data) {
#   # returns list of lists --
#   # each item of list contains input data for a mutation sample presence set
#   # original mutation indices from input_data are recorded in $mutation_indices
#   pres <- ifelse(input_data$y > 0, 1, 0)
#   pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
#   types <- sort(names(table(pat)), decreasing=TRUE)
#   if (length(types) == 1) {
#     type_indices <- list()
#     type_indices[[types]] <- seq_len(nrow(pres))
#   } else {
#     type_indices <- lapply(types, function(x) which(pat == x))
#     names(type_indices) <- types
#   }
#   warning('consider merge one-sample cluster with other clusters')
#   return(type_indices)
# }

#' separate copy number alterations by sample presence
# separateCNAsBySamplePresence <- function(input_data) {
#   # returns list of lists --
#   # each item of list contains input data for a mutation sample presence set
#   # original mutation indices from input_data are recorded in $mutation_indices
#   pres <- ifelse(input_data$tcn == 2, 0, 1)
#   pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
#   types <- sort(names(table(pat)), decreasing=TRUE)
#   if (length(types) == 1) {
#     type_indices <- list()
#     type_indices[[types]] <- seq_len(nrow(pres))
#   } else {
#     type_indices <- lapply(types, function(x) which(pat == x))
#     names(type_indices) <- types
#   }
#   return(type_indices)
# }