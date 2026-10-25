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
#' \code{noise_multiplier > max_noise}.
#'
#' Calibration (tree-delta, star graph, lambda = 0.01, 3-5 tumor samples and a normal sample; signed precision of
#' tumor-to-tumor pathway calls at padj < 0.05, i.e. pathway and direction both right; Roerink rows pooled over calls):
#' \tabular{lccc}{
#'   estimable fraction \tab < 0.8 \tab 0.8-0.95 \tab >= 0.95 \cr
#'   wellDR breast cancer, 8 patients (tree edges; subclones that differ little) \tab 0.06 \tab 0.18 \tab 0.49 \cr
#'   Roerink colorectal organoids, 3 patients (all 64 clone pairs; chance 0.18-0.27) \tab 0.41 \tab 0.67 \tab 0.76 \cr
#'   Roerink, calls with padj < 0.001 \tab 0.51 \tab 0.79 \tab 0.82 \cr
#' }
#' In both data sets precision rises with the estimable fraction. Where the true clone difference is small (wellDR),
#' calls on contrasts below 80\% estimable are almost all false; where clones differ strongly (Roerink), they are
#' still right more often than chance but less often than well-determined calls. This sets the default
#' \code{min_estimable = 0.95}: report calls on edges >= 0.95, report calls on 0.8-0.95 edges with a caveat, and treat
#' calls on edges < 0.8 as unsupported by the design. In the Roerink data (Poisson count noise), the noise multiplier
#' also mattered among edges >= 0.95 estimable: precision was 0.89, 0.78 and 0.51 for multipliers <= 20, 20-100 and
#' > 100 (in wellDR, without count noise, it did not change precision). This sets the default \code{max_noise = 100},
#' the cut-off used in the paper's tiers; \code{max_noise = Inf} reports the multiplier without using it for the flag.
#' For an edge that is not fully estimable, the noise multiplier covers only the estimable part of the contrast and is
#' therefore a lower bound.
#'
#' @param proportionFile,purityFile,treeFile As for \code{runDeconvolution()} (\code{treeFile}: the DNA tree whose
#'   edges are tested by \code{runGSEA()}).
#' @param normal_samples Names of RNA samples that are pure normal tissue (not in \code{proportionFile}).
#' @param min_estimable,max_noise Flag thresholds (defaults 0.95 and 100).
#' @return A data.frame with one row per edge.
#' @export
edgeResolvability <- function(proportionFile, treeFile, purityFile = NULL, normal_samples = character(0),
                              min_estimable = 0.95, max_noise = 100) {
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
