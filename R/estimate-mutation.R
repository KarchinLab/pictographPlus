#' estimate mutation cellular fraction
#' @export
#'
estimateMutation <- function(input_data, cna_est, cncf_est) {
  mcf_est <- input_data$y
  # for each box by mutation sample presence
  for (i in seq_len(nrow(input_data$y))) {
    # estimate mutation cellular fraction for each mutation
    # for (j in seq_len(length(sep_list$mutation[[i]]))) {
      mutation_index = i
      # if mutation overlaps cna
      if (any(input_data$overlap[mutation_index,]==1)) {
        
        # estimate mcf using CNA information
        cna_index = which(input_data$overlap[mutation_index,]==1)
        
        # known information
        m = 1 # multiplicity; can change
        y = input_data$y[mutation_index,] # alt read count; vector
        n = input_data$n[mutation_index,] # total read count; vector
        cncf = cncf_est[cna_index, ] # cncf; vector
        cna = cna_est[cna_index,] #cna: vector
        
        # estimate multiplicity such that all mcf <= 1
        m = estimate_multiplicity(y,n,cncf,cna)
        mcf =( y * (2 - 2 * cncf + cna * cncf) / n )- (m-1) * cncf
        mcf_est[mutation_index, ] = mcf
      } else {
        # estimate mcf assuming CNA == 2 and multiplicity == 1 in case mutation in cna neutral region
        mcf_est[mutation_index, ] = input_data$y[mutation_index,] * 2 / input_data$n[mutation_index,]
      }
    # }
  }
  return(mcf_est)
}

#' estimate mutation cellular fraction
#' export
#'
# estimateMutation1 <- function(input_data, sep_list, cna_est, cncf_est, m_est) {
#   mcf_est <- input_data$y
#   # for each box by mutation sample presence
#   for (i in seq_len(length(sep_list$mutation))) {
#     # estimate mutation cellular fraction for each mutation
#     for (j in seq_len(length(sep_list$mutation[[i]]))) {
#       mutation_index = sep_list$mutation[[i]][j]
#       # if mutation overlaps cna
#       if (any(input_data$overlap[mutation_index,]==1)) {
# 
#         # estimate mcf using CNA information
#         cna_index = which(input_data$overlap[mutation_index,]==1)
# 
#         # known information
#         m = 1 # multiplicity; can change
#         y = input_data$y[mutation_index,] # alt read count; vector
#         n = input_data$n[mutation_index,] # total read count; vector
#         cncf = cncf_est[mutation_index, ] # cncf; vector
#         cna = cna_est[mutation_index,] #cna: vector
#         m = m_est[mutation_index,]
#         # count number of cells arbitrarily
#         # cell_count = input_data$n[mutation_index,] / (2 - 2 * cncf_est[mutation_index, ] + cna_est[mutation_index,]* cncf_est[mutation_index, ])
# 
#         # case1: cna before mutation; multiplicity == 1
#         # case2: cna parallel of mutation; multiplicity == 1
#         # case3: cna after mutation; two proportions make up the mcf
# 
#         # check presence of mutation and cna
#         # tmp1 = ifelse(input_data$y[mutation_index,]>0,1,0)
#         # tmp2 = ifelse(input_data$tcn[cna_index,]==2,0,1)
#         # tmp3 = tmp1-tmp2
# 
#         # estimate multiplicity such that all mcf <= 1
#         # m = estimate_multiplicity(y,n,cncf,cna)
#         mcf =( y * (2 - 2 * cncf + cna * cncf) / n )- (m-1) * cncf
#         mcf_est[mutation_index, ] = mcf
#         # likelihood for case 1
# 
#         # likelihood for case 2
# 
#         # likelihood for case 3
# 
# 
#         # # which case is most likely
#         # if all(tmp3==0) {
#         #   # mutation and cna in the same box; all three cases possible
#         #   mcf = input_data$y[mutation_index,] * (2 - 2 * cncf_est[cna_index, ] + cna_est[cna_index,]* cncf_est[cna_index, ]) / input_data$n[mutation_index,] - (m-1) * cncf_est[cna_index, ]
#         # } else if (all(tmp3>=0)) {
#         #   # mutation in more samples than cna; case2 and case3 possible
#         # } else if (all(tmp3<=0)) {
#         #   # cna in more samples than mutation; case1 and case2 possible
#         # } else {
#         #   # mutation and cna samples not in order; case2 possible
#         # }
# 
# 
#       } else {
#         # estimate mcf assuming CNA == 2 and multiplicity == 1 in case mutation in cna neutral region
#         mcf_est[mutation_index, ] = input_data$y[mutation_index,] * 2 / input_data$n[mutation_index,]
#       }
#     }
#   }
#   return(mcf_est)
# }

#' estimate multiplicity for all mutations
#' @export
estimateMultiplicity1 <- function(input_data, cna_est, cncf_est) {
  m_est <- input_data$y
  # for each box by mutation sample presence
  for (i in seq_len(nrow(input_data$y))) {
    mutation_index = i
    # if mutation overlaps cna
    if (any(input_data$overlap[mutation_index,]==1)) {
      
      # estimate mcf using CNA information
      cna_index = which(input_data$overlap[mutation_index,]==1)
      
      # known information
      y = input_data$y[mutation_index,] # alt read count; vector
      n = input_data$n[mutation_index,] # total read count; vector
      cncf = cncf_est[mutation_index, ] # cncf; vector
      cna = cna_est[mutation_index,] #cna: vector

      # estimate multiplicity such that all mcf <= 1
      m = estimate_multiplicity(y,n,cncf,cna)
      m_est[mutation_index, ] = m
      } else {
        # estimate mcf assuming CNA == 2 and multiplicity == 1 in case mutation in cna neutral region
        m_est[mutation_index, ] = 1
      }
  }
  return(m_est)
}

#' estimate multiplicity for one gene such that all mcf <= 1
estimate_multiplicity <- function(y, n, cncf, cna, threshold=1) {
  warning("need to add copy number loss events")
  for (m in seq_len(cna[1])) {
    mcf = (y * (2 - 2 * cncf + cna * cncf) / n )- (m-1) * cncf
    if (all(mcf <= threshold)) {
      return(m)
    }
  }
  return(1)
}

# # calculate the likelihood for case1: cna before mutation; multiplicity == 1
# calcLikelihood_case1 <- function() {
# 
# }
# 
# # calculate the likelihood for case2: cna parallel of mutation; multiplicity == 1
# calcLikelihood_case2 <- function() {
# 
# }
# 
# # calculate the likelihood for case3: cna after mutation; two proportions make up the mcf
# calcLikelihood_case3 <- function() {
# 
# }

