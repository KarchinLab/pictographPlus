# ---- Edge resolvability (design diagnostic) ---------------------------
#
# An edge's pathway test compares two estimated clone profiles. How well the samples can separate those two clones
# depends only on the proportions: if the two clones rise and fall together across samples, their difference is
# poorly determined (or not determined at all, when there are fewer samples than clones) and the estimate is driven
# by regularization and noise. This diagnostic needs no expression data.

#' Per-edge resolvability of the clone contrast from the sample design
#'
#' Builds the mixing matrix exactly as \code{runDeconvolution()} does (tumor proportions x purity, non-tumor clone
#' "0" = 1 - purity, plus RNA-only normal samples as 100\% clone "0") and, for every tree edge, measures how well the
#' contrast child - parent is determined by it:
#' \itemize{
#'   \item \code{estimable_fraction}: the share (0-1) of the contrast that lies in the space the samples identify.
#'     Below 1, part of the edge's difference is not determined by the data at all and comes from regularization.
#'   \item \code{noise_multiplier}: \eqn{||a||^2 / 2}, where \eqn{a} weights each sample's expression in the
#'     least-squares (minimum-norm) estimate of the contrast. It is 1 when each clone has a pure sample of its own,
#'     and grows as the two clones' proportions become similar across samples.
#' }
#' An edge is flagged \code{resolvable = FALSE} when \code{estimable_fraction < min_estimable} or
#' \code{noise_multiplier > max_noise}. In the eight-patient wellDR benchmark, signed precision of tumor-to-tumor
#' pathway calls (tree-delta, star graph) was 0.06, 0.18 and 0.49 for edges < 80\%, 80-95\% and >= 95\% estimable,
#' which sets the default \code{min_estimable = 0.95}. The noise multiplier could not be calibrated there (the
#' benchmark mixtures carry no count noise), so by default it is reported but not used for the flag.
#'
#' @param proportionFile,purityFile,treeFile As for \code{runDeconvolution()} (\code{treeFile}: the DNA tree whose
#'   edges are tested by \code{runGSEA()}).
#' @param normal_samples Names of RNA samples that are pure normal tissue (not in \code{proportionFile}).
#' @param min_estimable,max_noise Flag thresholds (defaults 0.95 and \code{Inf}).
#' @return A data.frame with one row per edge.
#' @export
edgeResolvability <- function(proportionFile, treeFile, purityFile = NULL, normal_samples = character(0),
                              min_estimable = 0.95, max_noise = Inf) {
  prop <- as.matrix(read.csv(proportionFile, row.names = 1, check.names = FALSE))
  if (is.null(purityFile)) { pur <- rep(1, ncol(prop)) } else {
    pd <- read.csv(purityFile, check.names = FALSE); stopifnot(identical(colnames(pd), colnames(prop))); pur <- as.numeric(pd[1, ])
  }
  Pi <- rbind(`0` = 1 - pur, sweep(prop, 2, pur, "*"))
  for (s in normal_samples) Pi <- cbind(Pi, setNames(c(1, rep(0, nrow(prop))), rownames(Pi)))
  Pi <- t(Pi[order(as.integer(rownames(Pi))), , drop = FALSE])        # samples x clones, clone "0" first
  edges <- read_tree(treeFile); if (is.null(nrow(edges))) edges <- matrix(as.integer(edges), nrow = 1L)
  out <- .edge_resolvability(Pi, edges)
  out$resolvable <- out$estimable_fraction >= min_estimable & out$noise_multiplier <= max_noise
  out
}

# Pi: samples x clones (column k = clone k - 1, clone "0" first); edges: 0-based parent/child matrix.
.edge_resolvability <- function(Pi, edges) {
  sv <- svd(Pi); keep <- sv$d > max(dim(Pi)) * max(sv$d) * .Machine$double.eps
  M <- sv$v[, keep, drop = FALSE] %*% (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])   # pseudo-inverse: clones x samples
  Proj <- M %*% Pi                                                      # projector onto the identified clone space
  do.call(rbind, lapply(seq_len(nrow(edges)), function(k) {
    p <- edges[k, 1]; c <- edges[k, 2]; v <- numeric(ncol(Pi)); v[c + 1] <- 1; v[p + 1] <- -1
    a <- as.vector(v %*% M)
    data.frame(edge = sprintf("clone%d_%d", c, p), parent = p, child = c,
               estimable_fraction = round(sqrt(sum((Proj %*% v)^2) / sum(v^2)), 3),
               noise_multiplier = round(sum(a^2) / 2, 2))
  }))
}
