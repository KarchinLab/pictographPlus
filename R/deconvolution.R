#' Bulk RNA deconvolution using tumor evolution information
#'
#' @description
#' Deconvolves bulk RNA expression into clone-level profiles by integrating
#' tumor clonal tree structure and clone proportions. Seven model variants
#' are available; see the \code{model} parameter for details.
#'
#' @export
#' @import igraph
#'
#' @param rna_file CSV of raw gene counts; rows = genes, columns = samples.
#'   An optional normal sample column (not present in proportionFile) is
#'   automatically handled.
#' @param treeFile CSV tree file from \code{runPictograph} (or external tool).
#' @param proportionFile CSV subclone proportions from \code{runPictograph}.
#' @param outputDir Output directory for clonal_expression.csv.
#' @param normalize Normalize counts with DESeq2 size factors; default TRUE.
#' @param purityFile Optional purity CSV for tumor purity correction.
#' @param lambda Regularisation strength; default 0.01.
#' @param use_star_tree If TRUE (default), ignore \code{treeFile} and use a
#'   star topology (root directly connected to all clones). Recommended when
#'   the inferred tree topology does not improve performance. When FALSE,
#'   \code{treeFile} must be supplied.
#' @param model Deconvolution model. One of:
#'   \itemize{
#'     \item \code{"tree_delta"} (default) — group-L2 penalty on each edge's
#'       child-minus-parent difference, via ADMM. In the eight-patient wellDR-seq
#'       benchmark (lambda = 0.01, matched normal) it was the only model besides
#'       fused_ew to beat same-input NNLS on both edge-level pathway F1 and
#'       clone-specific (centered) recovery, with either the star or the true
#'       tree, and it was the most robust model in the operating-limit
#'       simulations. Slower than elastic_net on large panels (hours for six
#'       clones x 27k genes).
#'     \item \code{"elastic_net"} — L2 Laplacian + L1 fused LASSO via ADMM.
#'       Fast and reliably convergent; similar pathway F1 to tree_delta on the
#'       star graph but lower clone-specific recovery. The default before
#'       October 2026; requires lambda > 0.
#'     \item \code{"adaptive"} — iteratively reweighted Laplacian (IRLS).
#'       Highest MCC (0.248) in with-normal mode; use for low-FDR pathway calls.
#'     \item \code{"adaptive_v2"} — two-phase IRLS with unbiased initial
#'       weights. Best F1/sensitivity in tumor-only mode (no normal reference).
#'     \item \code{"plain"} — closed-form Laplacian smoothing. Best raw
#'       Pearson expression recovery.
#'     \item \code{"plain_debiased"} — plain Laplacian with shrinkage debiasing.
#'     \item \code{"fused_ew"} — element-wise fused LASSO via ADMM. Highest
#'       precision / MCC; conservative; requires lambda > 0.
#'   }
#' @param lambda_l2 Fixed L2 Laplacian weight for \code{"elastic_net"} model;
#'   default 0.01.
#' @param n_iter Number of IRLS iterations for adaptive models; default 5.
#' @param verbose Print solver convergence progress; default FALSE.
runDeconvolution <- function(rna_file,
                             treeFile       = NULL,
                             proportionFile,
                             outputDir,
                             normalize      = TRUE,
                             purityFile     = NULL,
                             lambda         = 0.01,
                             use_star_tree  = TRUE,
                             model          = "tree_delta",
                             lambda_l2      = 0.01,
                             n_iter         = 5,
                             verbose        = FALSE) {

  valid_models <- c("plain", "adaptive", "adaptive_v2", "plain_debiased",
                    "fused_ew", "elastic_net", "tree_delta")
  if (!model %in% valid_models) {
    stop(sprintf("Unknown model '%s'. Choose from: %s",
                 model, paste(valid_models, collapse = ", ")))
  }
  if (model %in% c("fused_ew", "elastic_net") && lambda <= 0) {
    stop(sprintf("Model '%s' requires lambda > 0.", model))
  }

  rnaData <- read.csv(rna_file, row.names = 1)
  if (normalize) {
    rnaData <- normalize_RNA(rnaData)
  }

  propData <- read.csv(proportionFile, row.names = 1)

  if (is.null(purityFile)) {
    proportionDF <- propData
    proportionDF <- rbind(proportionDF, 0)
  } else {
    purityData <- read.csv(purityFile)
    # purityData[1, ] is matched to propData's columns by POSITION (as.numeric()
    # drops names); assert the order actually agrees instead of silently
    # multiplying the wrong purity onto the wrong sample (see case_config.R).
    stopifnot(identical(colnames(purityData), colnames(propData)))
    proportionDF <- t(t(propData) * as.numeric(purityData[1, ]))
    proportionDF <- rbind(proportionDF, 1 - purityData)
  }

  rownames(proportionDF)[nrow(proportionDF)] <- "0"

  # Add normal sample column if present in RNA but absent from proportions
  lesion <- setdiff(colnames(rnaData), colnames(proportionDF))
  if (length(lesion) > 0) {
    for (s in lesion) {
      proportionDF[[s]] <- 0
      proportionDF[nrow(proportionDF), s] <- 1
    }
  }

  proportionDF <- round(proportionDF, 4)
  proportionDF <- proportionDF[, colnames(rnaData), drop = FALSE]
  # Numeric sort: rownames are clone-id strings ("0" = normal), and a
  # lexicographic sort misorders them once there are >=10 tumor clones
  # (e.g. "10" before "2"), which silently mis-indexes edges built from
  # treeFile for use_star_tree = FALSE.
  proportionDF <- proportionDF[order(as.integer(rownames(proportionDF))), , drop = FALSE]

  # Build 0-based edge list (star tree or from file)
  if (use_star_tree) {
    leaf_ids <- sort(as.integer(setdiff(rownames(proportionDF), "0")))
    n_leaves <- length(leaf_ids)
    edges    <- matrix(c(rep(0L, n_leaves), leaf_ids), ncol = 2L)
    storage.mode(edges) <- "integer"
    star_df  <- data.frame(edge_id = paste0("root->", leaf_ids),
                            parent  = "root",
                            child   = as.character(leaf_ids))
    write.csv(star_df, file.path(outputDir, "star_tree.csv"),
              row.names = FALSE, quote = FALSE)
  } else {
    if (is.null(treeFile)) stop("treeFile must be supplied when use_star_tree = FALSE")
    edges <- read_tree(treeFile)
    if (is.null(nrow(edges))) {
      edges <- matrix(as.integer(edges), nrow = 1L)
    }
  }

  Y  <- as.matrix(t(rnaData))       # n_samples x n_genes
  Pi <- as.matrix(t(proportionDF))  # n_samples x K_clones

  # Dispatch to chosen model
  fit <- switch(model,
    plain = fit_plain_laplacian(Y, Pi, edges, lambda = lambda, verbose = verbose),
    adaptive = fit_adaptive_laplacian(Y, Pi, edges, lambda = lambda,
                                      n_iter = n_iter, verbose = verbose),
    adaptive_v2 = fit_adaptive_v2(Y, Pi, edges, lambda = lambda,
                                   n_iter = n_iter, verbose = verbose),
    plain_debiased = {
      fp <- fit_plain_laplacian(Y, Pi, edges, lambda = lambda, verbose = verbose)
      list(X = debias_laplacian(fp$X, Pi, edges, lambda = lambda))
    },
    fused_ew = fit_elementwise_fused_lasso_admm(Y, Pi, edges, lambda = lambda,
                                                 verbose = verbose),
    elastic_net = fit_elastic_net_tree(Y, Pi, edges, lambda1 = lambda_l2,
                                        lambda2 = lambda, verbose = verbose),
    tree_delta  = fit_tree_delta_admm(Y, Pi, edges, lambda = lambda,
                                       verbose = verbose)
  )

  X_optimal <- fit$X
  rownames(X_optimal) <- colnames(Pi)   # clone IDs as row names
  colnames(X_optimal) <- colnames(Y)    # gene names as column names

  write.csv(X_optimal, file = file.path(outputDir, "clonal_expression.csv"))
  # Solver record (ADMM models report convergence; closed-form models leave these NA).
  fit_info <- data.frame(model = model, lambda = lambda, use_star_tree = use_star_tree,
                         n_samples = nrow(Y), n_clones = ncol(Pi), n_genes = ncol(Y),
                         converged = if (is.null(fit$converged)) NA else fit$converged,
                         exit_reason = if (is.null(fit$exit_reason)) NA else fit$exit_reason,
                         iterations = if (is.null(fit$iterations)) NA else fit$iterations,
                         x_kkt_residual = if (is.null(fit$x_kkt_residual)) NA else fit$x_kkt_residual)
  write.csv(fit_info, file = file.path(outputDir, "deconvolution_fit_info.csv"), row.names = FALSE)
  if (isFALSE(fit$converged))
    warning(sprintf("%s did not converge within max_iter; see deconvolution_fit_info.csv", model))
  return(X_optimal)
}


# ---- DESeq2 normalisation -------------------------------------------

#' @import DESeq2
normalize_RNA <- function(rnaData) {
  counts  <- as.matrix(rnaData)
  colData <- data.frame(row.names = colnames(counts))
  dds     <- DESeqDataSetFromMatrix(countData = counts, colData = colData,
                                    design = ~ 1)
  dds     <- estimateSizeFactors(dds)
  as.data.frame(counts(dds, normalized = TRUE))
}


# ---- Tree I/O -------------------------------------------------------

#' Read tree edge list (0-based integer matrix)
#' @export
read_tree <- function(treeFile) {
  edges <- list()
  con   <- file(treeFile, "r")
  on.exit(close(con))
  readLines(con, n = 1)   # skip header

  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    parts  <- strsplit(line, ",")[[1]]
    parent <- trimws(parts[2])
    child  <- trimws(parts[3])
    parent_int <- if (parent == "root") 0L else as.integer(parent)
    child_int  <- if (child  == "root") 0L else as.integer(child)
    edges <- append(edges, list(c(parent_int, child_int)))
  }

  edge_mat <- as.matrix(do.call(rbind, edges))
  storage.mode(edge_mat) <- "integer"
  edge_mat
}


# ---- Laplacian helpers ----------------------------------------------

build_laplacian <- function(edges, K, weights = NULL, normalize = "spectral") {
  if (is.null(weights)) weights <- rep(1.0, nrow(edges))
  L <- matrix(0.0, K, K)
  for (e in seq_len(nrow(edges))) {
    i <- edges[e, 1L] + 1L
    j <- edges[e, 2L] + 1L
    w <- weights[e]
    L[i, i] <- L[i, i] + w
    L[j, j] <- L[j, j] + w
    L[i, j] <- L[i, j] - w
    L[j, i] <- L[j, i] - w
  }
  if (normalize == "spectral") {
    d          <- pmax(diag(L), 1e-10)
    D_inv_sqrt <- diag(1.0 / sqrt(d))
    L          <- D_inv_sqrt %*% L %*% D_inv_sqrt
  }
  L
}

build_incidence <- function(edges, K) {
  n_edges <- nrow(edges)
  D <- matrix(0.0, n_edges, K)
  for (e in seq_len(n_edges)) {
    D[e, edges[e, 1L] + 1L] <- -1.0
    D[e, edges[e, 2L] + 1L] <-  1.0
  }
  D
}


# ---- Model 1: Plain Laplacian (closed-form) -------------------------

fit_plain_laplacian <- function(Y, Pi, edges, lambda = 0.05, ridge = 1e-8,
                                 normalize = "spectral", verbose = FALSE) {
  K <- max(edges) + 1L
  L <- build_laplacian(edges, K, normalize = normalize)
  A <- t(Pi) %*% Pi + lambda * L + ridge * diag(K)
  B <- t(Pi) %*% Y
  X <- tryCatch(solve(A, B), error = function(e) qr.solve(A, B))
  X <- pmax(X, 0)
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"))
}


# ---- Model 2: Adaptive Laplacian (IRLS) -----------------------------

fit_adaptive_laplacian <- function(Y, Pi, edges, lambda = 0.05, ridge = 1e-8,
                                    normalize = "spectral", n_iter = 5,
                                    eps = 1e-6, verbose = FALSE) {
  K       <- max(edges) + 1L
  weights <- rep(1.0, nrow(edges))
  X       <- NULL

  for (iter in seq_len(n_iter)) {
    L <- build_laplacian(edges, K, weights = weights, normalize = normalize)
    A <- t(Pi) %*% Pi + lambda * L + ridge * diag(K)
    B <- t(Pi) %*% Y
    X <- tryCatch(solve(A, B), error = function(e) qr.solve(A, B))
    X <- pmax(X, 0)
    for (e in seq_len(nrow(edges))) {
      i          <- edges[e, 1L] + 1L
      j          <- edges[e, 2L] + 1L
      diff_norm  <- sqrt(sum((X[i, ] - X[j, ])^2))
      weights[e] <- 1.0 / (diff_norm + eps)
    }
    if (verbose) cat(sprintf("  Adaptive IRLS iter %d\n", iter))
  }
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"))
}


# ---- Model 3: Adaptive v2 (two-phase IRLS) --------------------------

fit_adaptive_v2 <- function(Y, Pi, edges, lambda = 0.05, ridge = 1e-8,
                              normalize = "spectral", n_iter = 5,
                              eps = 1e-4, verbose = FALSE) {
  K       <- max(edges) + 1L
  n_edges <- nrow(edges)

  # Phase 1: unregularised solve for unbiased initial weights
  A0     <- t(Pi) %*% Pi + ridge * diag(K)
  X_init <- tryCatch(solve(A0, t(Pi) %*% Y),
                     error = function(e) qr.solve(A0, t(Pi) %*% Y))
  X_init <- pmax(X_init, 0)

  # Phase 2: initial edge weights from unbiased differences
  weights <- numeric(n_edges)
  for (e in seq_len(n_edges)) {
    i          <- edges[e, 1L] + 1L
    j          <- edges[e, 2L] + 1L
    diff_norm  <- sqrt(sum((X_init[i, ] - X_init[j, ])^2))
    weights[e] <- 1.0 / (diff_norm + eps)
  }

  # Phase 3: IRLS with informed weights
  X <- X_init
  for (iter in seq_len(n_iter)) {
    L <- build_laplacian(edges, K, weights = weights, normalize = normalize)
    A <- t(Pi) %*% Pi + lambda * L + ridge * diag(K)
    X <- tryCatch(solve(A, t(Pi) %*% Y),
                  error = function(e) qr.solve(A, t(Pi) %*% Y))
    X <- pmax(X, 0)
    for (e in seq_len(n_edges)) {
      i          <- edges[e, 1L] + 1L
      j          <- edges[e, 2L] + 1L
      diff_norm  <- sqrt(sum((X[i, ] - X[j, ])^2))
      weights[e] <- 1.0 / (diff_norm + eps)
    }
    if (verbose) cat(sprintf("  Adaptive-v2 IRLS iter %d\n", iter))
  }
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"))
}


# ---- Model 4: Laplacian shrinkage debiasing -------------------------
# X_debiased = X_pen + lambda * (Pi'Pi + rI)^{-1} * L * X_pen

debias_laplacian <- function(X_pen, Pi, edges, lambda, ridge = 1e-8,
                               normalize = "spectral") {
  K   <- max(edges) + 1L
  L   <- build_laplacian(edges, K, normalize = normalize)
  PtP <- t(Pi) %*% Pi
  sv  <- svd(PtP)

  tol  <- max(sv$d) * 0.01
  keep <- sv$d > tol
  if (sum(keep) == 0) return(X_pen)

  d_inv <- numeric(length(sv$d))
  d_inv[keep] <- 1.0 / (sv$d[keep] + ridge)
  PtP_pinv    <- sv$v %*% diag(d_inv) %*% t(sv$u)

  pmax(X_pen + lambda * PtP_pinv %*% L %*% X_pen, 0)
}


# ---- Shared X-subproblem solver (validated_admm_solver_v3.R) -------
#
# Each ADMM model below solves the same kind of X-subproblem each outer
# iteration: min_{X>=0} X'HX - 2 B'X for a fixed PSD H. The historical
# implementation computed the unconstrained minimizer H^{-1}B and clipped
# negative entries with pmax(.,0) -- not the constrained minimizer unless H
# is diagonal, and on real panels this stopped at an objective >2x worse
# than the correct solve (analysis/AUDIT_2026-09-15_TECHNICAL_VALIDATION.md).
# `.projected_qp` instead solves it by projected FISTA, with an active-set
# KKT polish and a recorded KKT residual so callers can confirm convergence.
.projected_qp <- function(H, B, X0 = NULL, tol = 1e-8, max_iter = 5000L, kkt_tol = 1e-5) {
  K <- nrow(H)
  X <- if (is.null(X0)) matrix(0, K, ncol(B)) else pmax(X0, 0)
  Yk <- X
  tk <- 1
  lipschitz <- max(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(lipschitz) || lipschitz <= 0) stop("X-subproblem Hessian is not positive definite")
  converged <- FALSE
  rel_change <- Inf
  for (it in seq_len(max_iter)) {
    Xnew <- pmax(Yk - (H %*% Yk - B) / lipschitz, 0)
    rel_change <- norm(Xnew - X, "F") / (norm(X, "F") + 1e-12)
    if (rel_change < tol) {
      grad_new <- H %*% Xnew - B
      kkt_new <- max(abs(ifelse(Xnew > 1e-10, grad_new, pmin(grad_new, 0))))
      if (kkt_new > kkt_tol) {
        active <- Xnew > 1e-10
        groups <- split(seq_len(ncol(B)), apply(active, 2, paste0, collapse = ""))
        for (cols in groups) {
          free <- which(active[, cols[1]])
          candidate <- matrix(0, K, length(cols))
          if (length(free)) {
            sol <- tryCatch(solve(H[free, free, drop = FALSE], B[free, cols, drop = FALSE]),
                            error = function(e) NULL)
            if (is.null(sol)) next
            candidate[free, ] <- sol
          }
          cg <- H %*% candidate - B[, cols, drop = FALSE]
          ck <- apply(abs(ifelse(candidate > 1e-10, cg, pmin(cg, 0))), 2, max)
          accept <- colSums(candidate < 0) == 0 & is.finite(ck) & ck <= kkt_tol
          if (any(accept)) Xnew[, cols[accept]] <- candidate[, accept, drop = FALSE]
        }
        grad_new <- H %*% Xnew - B
        kkt_new <- max(abs(ifelse(Xnew > 1e-10, grad_new, pmin(grad_new, 0))))
      }
      if (kkt_new <= kkt_tol) { X <- Xnew; converged <- TRUE; break }
    }
    tnew <- (1 + sqrt(1 + 4 * tk^2)) / 2
    Yk <- Xnew + ((tk - 1) / tnew) * (Xnew - X)
    X <- Xnew
    tk <- tnew
  }
  grad <- H %*% X - B
  kkt <- max(abs(ifelse(X > 1e-10, grad, pmin(grad, 0))))
  list(X = X, iterations = it, converged = converged,
       rel_change = rel_change, kkt_residual = kkt)
}

# `.exact_qp` solves the same X-subproblem exactly, column by column, by block
# principal pivoting NNLS (Kim & Park 2011, SIAM J Sci Comput 33:3261), warm-
# started from the previous ADMM iteration's passive set `F0` (logical K x n).
# Columns sharing a passive set share one Cholesky solve. It is what the ADMM
# fits below use: same iterates as `.projected_qp`, 10-25x faster per
# iteration on the case studies, and exact where FISTA stalls on ill-conditioned
# H (analysis/solver_speed_2026-10-08/). Any column that does not settle within
# `max_it` pivots falls back to `.projected_qp`.
.exact_qp <- function(H, B, F0, max_it = 50L, kkt_tol = 1e-5,
                      x_tol = 1e-8, x_max_iter = 5000L) {
  K <- nrow(H); n <- ncol(B)
  w <- 2^(seq_len(K) - 1)
  F <- F0
  X <- matrix(0, K, n)
  solve_cols <- function(cols) {
    code <- colSums(F[, cols, drop = FALSE] * w)
    for (cc in unique(code)) {
      cj <- cols[code == cc]; free <- which(F[, cj[1]])
      X[, cj] <<- 0
      if (length(free)) {
        R <- tryCatch(chol(H[free, free, drop = FALSE]), error = function(e) NULL)
        X[free, cj] <<- if (is.null(R)) qr.solve(H[free, free, drop = FALSE], B[free, cj, drop = FALSE])
                        else backsolve(R, forwardsolve(t(R), B[free, cj, drop = FALSE]))
      }
    }
  }
  solve_cols(seq_len(n))
  # gradient tolerance for "active and optimal"; kept well below kkt_tol so the
  # outer KKT gate can pass on large-count data (max|B| ~ 1e8)
  ytol <- min(1e-12 * max(1, max(abs(B))), 0.01 * kkt_tol)
  alpha <- rep(3L, n); beta <- rep(K + 1L, n)
  todo <- seq_len(n)
  for (it in seq_len(max_it)) {
    Xt <- X[, todo, drop = FALSE]; Ft <- F[, todo, drop = FALSE]
    G  <- H %*% Xt - B[, todo, drop = FALSE]
    inf <- (Ft & Xt < 0) | (!Ft & G < -ytol)
    ninf <- colSums(inf)
    todo_new <- todo[ninf > 0]
    if (!length(todo_new)) { todo <- integer(0); break }
    inf <- inf[, ninf > 0, drop = FALSE]; ninf <- ninf[ninf > 0]
    full <- ninf < beta[todo_new] | alpha[todo_new] >= 1L
    dec  <- ninf < beta[todo_new]
    beta[todo_new[dec]]  <- ninf[dec]; alpha[todo_new[dec]] <- 3L
    nd <- !dec & full; alpha[todo_new[nd]] <- alpha[todo_new[nd]] - 1L
    flip <- inf
    if (any(!full)) {   # backup rule: flip only the largest infeasible index
      for (j in which(!full)) { r <- max(which(inf[, j])); flip[, j] <- FALSE; flip[r, j] <- TRUE }
    }
    Fsub <- F[, todo_new, drop = FALSE]; Fsub[flip] <- !Fsub[flip]; F[, todo_new] <- Fsub
    todo <- todo_new
    solve_cols(todo)
  }
  if (length(todo)) {
    fb <- .projected_qp(H, B[, todo, drop = FALSE], pmax(X[, todo, drop = FALSE], 0),
                        tol = x_tol, max_iter = x_max_iter, kkt_tol = kkt_tol)
    X[, todo] <- fb$X
  }
  X <- pmax(X, 0)       # clears -0 / 1e-17 round-off on the passive set
  grad <- H %*% X - B
  list(X = X, F = X > 0, n_fallback = length(todo),
       kkt_residual = max(abs(ifelse(X > 1e-10, grad, pmin(grad, 0)))))
}


# ---- Model 5: Element-wise fused LASSO via ADMM --------------------

fit_elementwise_fused_lasso_admm <- function(Y, Pi, edges, lambda = 0.01,
                                              rho = 1.0, max_iter = 60000L,
                                              tol = 1e-6, ridge = 1e-8,
                                              adaptive_rho = TRUE,
                                              rho_mu = 10.0, rho_tau = 2.0,
                                              x_tol = 1e-8, x_max_iter = 5000L,
                                              kkt_tol = 1e-5, x_stall_tol = 0,
                                              verbose = FALSE) {
  # `adaptive_rho` turns on Boyd 2011 sec 3.4.1 residual balancing: every 25
  # iterations from iteration 100, if the primal / dual residuals have drifted
  # > `rho_mu`x apart, rho is nudged by `rho_tau` and the scaled dual u rescaled
  # with it. rho is not bounded and adaptation never freezes: the earlier
  # [rho0/64, rho0*64] bounds, 50% freeze and 40-change cap kept under-determined
  # tree_delta fits from converging (analysis/solver_speed_2026-10-08/). Without it the
  # fixed-rho iteration stalls with a flat dual residual on under-determined Pi
  # and hits max_iter without tripping `tol` (see analysis/case_study_model_switch/,
  # test c), which invents spurious distal-edge GSEA signal. Updating rho only
  # periodically (not every iteration) lets ADMM equilibrate between changes.
  # The X-update is the exact non-negative solve (`.exact_qp` above);
  # `tol` sets both the absolute and relative Boyd stopping tolerance for the
  # outer ADMM primal/dual residuals, which must also clear the X-subproblem's
  # own KKT gate (`kkt_tol`) before an iteration counts as converged.
  # `adaptive_rho = FALSE` + `max_iter = 5000` + `tol = 1e-3` reproduces the
  # pre-2026-09 iteration budget, not the historical fits (those clipped
  # instead of solving the X-subproblem; see the header note above).
  K       <- max(edges) + 1L
  n_genes <- ncol(Y)
  n_edges <- nrow(edges)
  D       <- build_incidence(edges, K)
  Dt      <- t(D)
  DtD     <- Dt %*% D
  PtP2    <- 2.0 * t(Pi) %*% Pi
  PtY2    <- 2.0 * t(Pi) %*% Y

  make_H <- function(r) PtP2 + r * DtD + ridge * diag(K)
  H <- make_H(rho)

  X <- matrix(0.0, K, n_genes)
  Z <- matrix(0.0, n_edges, n_genes)
  u <- matrix(0.0, n_edges, n_genes)
  iters <- 0L
  converged <- FALSE
  exit_reason <- NA_character_  # "kkt_gate" (primal/dual/KKT all in tolerance) or "x_stall" (X
                                 # plateaued but the KKT gate was not necessarily met) or NA (hit
                                 # max_iter) -- readiness-review finding 2, 2026-09-21: callers must
                                 # be able to tell which exit fired, since only "kkt_gate" is the
                                 # strict convergence certificate other gated pipelines use.
  p_res <- d_res <- NA_real_
  n_rho_upd  <- 0L
  Fset <- matrix(TRUE, K, n_genes)   # warm-start passive set for .exact_qp
  X_ref <- X; chk_every <- 200L   # objective-plateau early stop (x_stall_tol <= 0, the default, disables it: only the primal/dual/KKT gate stops the fit)
  xinfo <- NULL

  for (iter in seq_len(max_iter)) {
    iters  <- iter
    xinfo  <- .exact_qp(H, PtY2 + rho * Dt %*% (Z - u), Fset, kkt_tol = kkt_tol,
                        x_tol = x_tol, x_max_iter = x_max_iter)
    Fset <- xinfo$F
    X_new  <- xinfo$X
    V      <- D %*% X_new + u
    Z_new  <- sign(V) * pmax(abs(V) - lambda / rho, 0)
    resid  <- D %*% X_new - Z_new
    u      <- u + resid

    DX_n   <- norm(D %*% X_new, "F"); Zn <- norm(Z_new, "F")
    p_res  <- norm(resid, "F")
    d_res  <- norm(rho * Dt %*% (Z_new - Z), "F")
    eps_pri  <- sqrt(length(resid)) * tol + tol * max(DX_n, Zn)
    eps_dual <- sqrt(length(X_new)) * tol + tol * norm(rho * Dt %*% u, "F")
    Z <- Z_new; X <- X_new

    if (verbose && iter %% 100 == 0)
      cat(sprintf("  fused_ew ADMM iter %d: r=%.2e/%.2e s=%.2e/%.2e xKKT=%.2e rho=%.3g\n",
                  iter, p_res, eps_pri, d_res, eps_dual, xinfo$kkt_residual, rho))
    if (p_res <= eps_pri && d_res <= eps_dual && xinfo$kkt_residual <= kkt_tol) {
      converged <- TRUE; exit_reason <- "kkt_gate"; break
    }
    # once the primal constraint holds, the dual can crawl in the Pi null space
    # for a very long time while X is already stationary -- stop then.
    if (x_stall_tol > 0 && iter %% chk_every == 0L) {
      if (norm(X - X_ref, "F") / (norm(X_ref, "F") + 1e-12) < x_stall_tol) {
        converged <- TRUE; exit_reason <- "x_stall"; break
      }
      X_ref <- X
    }

    if (adaptive_rho && iter %% 25L == 0L && iter >= 100L) {
      new_rho <- if (p_res > rho_mu * d_res) rho * rho_tau
                 else if (d_res > rho_mu * p_res) rho / rho_tau
                 else rho
      if (new_rho != rho) {
        u <- u * (rho / new_rho); rho <- new_rho
        H <- make_H(rho)
        n_rho_upd <- n_rho_upd + 1L
      }
    }
  }
  if (verbose && !converged)
    cat(sprintf("  fused_ew ADMM: hit max_iter=%d (r=%.2e s=%.2e xKKT=%.2e)\n",
                max_iter, p_res, d_res, xinfo$kkt_residual))
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"),
       iterations = iters, rho_final = rho, converged = converged,
       exit_reason = exit_reason, x_kkt_residual = xinfo$kkt_residual,
       primal_residual = p_res, dual_residual = d_res,
       eps_primal = eps_pri, eps_dual = eps_dual)
}


# ---- Model 6: Elastic net (L2 + L1) via ADMM -----------------------

fit_elastic_net_tree <- function(Y, Pi, edges, lambda1 = 0.01, lambda2 = 0.01,
                                  ridge = 1e-8, normalize = "spectral",
                                  max_iter = 60000L, tol = 1e-6,
                                  adaptive_rho = TRUE,
                                  rho_mu = 10.0, rho_tau = 2.0,
                                  x_tol = 1e-8, x_max_iter = 5000L,
                                  kkt_tol = 1e-5, x_stall_tol = 0,
                                  verbose = FALSE) {
  # See fit_elementwise_fused_lasso_admm for `adaptive_rho`, the X-subproblem
  # solve (`.exact_qp`), and what `tol`/`kkt_tol` gate. The L2 Laplacian
  # term here already conditions the ADMM system, so this solver converged even
  # under fixed rho; residual balancing mainly speeds it up.
  K       <- max(edges) + 1L
  n_genes <- ncol(Y)
  n_edges <- nrow(edges)
  L       <- build_laplacian(edges, K, normalize = normalize)
  D       <- build_incidence(edges, K)
  Dt      <- t(D)
  DtD     <- Dt %*% D
  PtP2    <- 2.0 * t(Pi) %*% Pi
  PtY2    <- 2.0 * t(Pi) %*% Y

  rho   <- 1.0
  make_H <- function(r) PtP2 + lambda1 * L + r * DtD + ridge * diag(K)
  H <- make_H(rho)

  X <- matrix(0.0, K, n_genes)
  Z <- matrix(0.0, n_edges, n_genes)
  u <- matrix(0.0, n_edges, n_genes)
  iters <- 0L
  converged <- FALSE
  exit_reason <- NA_character_  # "kkt_gate" (primal/dual/KKT all in tolerance) or "x_stall" (X
                                 # plateaued but the KKT gate was not necessarily met) or NA (hit
                                 # max_iter) -- readiness-review finding 2, 2026-09-21: callers must
                                 # be able to tell which exit fired, since only "kkt_gate" is the
                                 # strict convergence certificate other gated pipelines use.
  p_res <- d_res <- NA_real_
  n_rho_upd  <- 0L
  Fset <- matrix(TRUE, K, n_genes)   # warm-start passive set for .exact_qp
  X_ref <- X; chk_every <- 200L
  xinfo <- NULL

  for (iter in seq_len(max_iter)) {
    iters  <- iter
    xinfo  <- .exact_qp(H, PtY2 + rho * Dt %*% (Z - u), Fset, kkt_tol = kkt_tol,
                        x_tol = x_tol, x_max_iter = x_max_iter)
    Fset <- xinfo$F
    X_new  <- xinfo$X
    V      <- D %*% X_new + u
    Z_new  <- sign(V) * pmax(abs(V) - lambda2 / rho, 0)
    resid  <- D %*% X_new - Z_new
    u      <- u + resid

    DX_n  <- norm(D %*% X_new, "F"); Zn <- norm(Z_new, "F")
    p_res <- norm(resid, "F")
    d_res <- norm(rho * Dt %*% (Z_new - Z), "F")
    eps_pri  <- sqrt(length(resid)) * tol + tol * max(DX_n, Zn)
    eps_dual <- sqrt(length(X_new)) * tol + tol * norm(rho * Dt %*% u, "F")
    Z <- Z_new; X <- X_new

    if (verbose && iter %% 100 == 0)
      cat(sprintf("  elastic_net ADMM iter %d: r=%.2e/%.2e s=%.2e/%.2e xKKT=%.2e rho=%.3g\n",
                  iter, p_res, eps_pri, d_res, eps_dual, xinfo$kkt_residual, rho))
    if (p_res <= eps_pri && d_res <= eps_dual && xinfo$kkt_residual <= kkt_tol) {
      converged <- TRUE; exit_reason <- "kkt_gate"; break
    }
    if (x_stall_tol > 0 && iter %% chk_every == 0L) {
      if (norm(X - X_ref, "F") / (norm(X_ref, "F") + 1e-12) < x_stall_tol) {
        converged <- TRUE; exit_reason <- "x_stall"; break
      }
      X_ref <- X
    }

    if (adaptive_rho && iter %% 25L == 0L && iter >= 100L) {
      new_rho <- if (p_res > rho_mu * d_res) rho * rho_tau
                 else if (d_res > rho_mu * p_res) rho / rho_tau
                 else rho
      if (new_rho != rho) {
        u <- u * (rho / new_rho); rho <- new_rho
        H <- make_H(rho)
        n_rho_upd <- n_rho_upd + 1L
      }
    }
  }
  if (verbose && !converged)
    cat(sprintf("  elastic_net ADMM: hit max_iter=%d (r=%.2e s=%.2e xKKT=%.2e)\n",
                max_iter, p_res, d_res, xinfo$kkt_residual))
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"),
       iterations = iters, rho_final = rho, converged = converged,
       exit_reason = exit_reason, x_kkt_residual = xinfo$kkt_residual,
       primal_residual = p_res, dual_residual = d_res,
       eps_primal = eps_pri, eps_dual = eps_dual)
}


# ---- Model 7: Tree-delta parameterisation via ADMM -----------------
# X_c = mu + sum_{e on path(root->c)} delta_e
# Penalty: lambda * sum_e ||delta_e||_2  (group L2, promotes sparsity over edges)

build_path_matrix <- function(edges, K) {
  n_edges <- nrow(edges)
  stopifnot(n_edges == K - 1L)

  children <- vector("list", K)
  for (e in seq_len(n_edges)) {
    p <- edges[e, 1L] + 1L
    ch <- edges[e, 2L] + 1L
    children[[p]] <- c(children[[p]], ch)
  }

  all_ch  <- edges[, 2L] + 1L
  root_1b <- setdiff(edges[, 1L] + 1L, all_ch)
  if (length(root_1b) != 1L) root_1b <- 1L

  T_mat <- matrix(0.0, K, K)

  dfs <- function(node, active_edges) {
    T_mat[node, 1L] <<- 1.0
    for (ei in active_edges) T_mat[node, 1L + ei] <<- 1.0
    for (ch in children[[node]]) {
      ei <- which(edges[, 1L] == (node - 1L) & edges[, 2L] == (ch - 1L))
      dfs(ch, c(active_edges, ei))
    }
  }
  dfs(root_1b, integer(0L))
  T_mat
}

fit_tree_delta_admm <- function(Y, Pi, edges, lambda = 0.05,
                                 rho = 1.0, max_iter = 60000L, tol = 1e-6,
                                 ridge = 1e-8, adaptive_rho = TRUE,
                                 rho_mu = 10.0, rho_tau = 2.0,
                                 x_tol = 1e-8, x_max_iter = 5000L,
                                 kkt_tol = 1e-5, x_stall_tol = 0, verbose = FALSE) {
  # See fit_elementwise_fused_lasso_admm for `adaptive_rho` and the X-subproblem
  # solve (`.exact_qp`). On a rooted tree, an edge's row of D %*% X is
  # exactly delta_e (child expression minus parent expression), so the group-L2
  # penalty is solved directly in X coordinates via the shared ADMM/QP machinery
  # instead of the historical Delta/T_mat reparameterisation, whose closed-form
  # (lambda = 0) branch and post-hoc `pmax(T_mat %*% Delta, 0)` projection never
  # solved the constrained problem (on a real P8 panel: objective 27604.9
  # historical vs 30.5 corrected; see
  # analysis/case_study_model_switch/results/real_benchmark_admm_validation.csv).
  # `Delta`/`T_mat` are still returned (some callers persist `Delta_opt.rds`),
  # recovered from the converged X: `T_mat` is invertible for any rooted tree,
  # so `Delta = solve(T_mat, X)` recovers root expression (row 1) and each
  # edge's delta (row 1+e, matching `edges`' row order) exactly.
  K       <- max(edges) + 1L
  n_genes <- ncol(Y)
  n_edges <- nrow(edges)
  T_mat   <- build_path_matrix(edges, K)
  D       <- build_incidence(edges, K)
  Dt      <- t(D)
  DtD     <- Dt %*% D
  PtP2    <- 2.0 * t(Pi) %*% Pi
  PtY2    <- 2.0 * t(Pi) %*% Y

  make_H <- function(r) PtP2 + r * DtD + ridge * diag(K)
  H <- make_H(rho)

  X <- matrix(0.0, K, n_genes)
  Z <- matrix(0.0, n_edges, n_genes)
  u <- matrix(0.0, n_edges, n_genes)
  iters <- 0L
  converged <- FALSE
  exit_reason <- NA_character_  # "kkt_gate" (primal/dual/KKT all in tolerance) or "x_stall" (X
                                 # plateaued but the KKT gate was not necessarily met) or NA (hit
                                 # max_iter) -- readiness-review finding 2, 2026-09-21: callers must
                                 # be able to tell which exit fired, since only "kkt_gate" is the
                                 # strict convergence certificate other gated pipelines use.
  p_res <- d_res <- NA_real_
  n_rho_upd  <- 0L
  Fset <- matrix(TRUE, K, n_genes)   # warm-start passive set for .exact_qp
  X_ref <- X; chk_every <- 200L
  xinfo <- NULL

  for (iter in seq_len(max_iter)) {
    iters <- iter
    xinfo <- .exact_qp(H, PtY2 + rho * Dt %*% (Z - u), Fset, kkt_tol = kkt_tol,
                       x_tol = x_tol, x_max_iter = x_max_iter)
    Fset <- xinfo$F
    X_new <- xinfo$X

    V     <- D %*% X_new + u
    Z_new <- matrix(0.0, n_edges, n_genes)
    vn    <- sqrt(rowSums(V^2))
    active <- vn > lambda / rho
    if (any(active))
      Z_new[active, ] <- V[active, , drop = FALSE] * (1 - lambda / rho / vn[active])

    resid <- D %*% X_new - Z_new
    u     <- u + resid

    DX_n  <- norm(D %*% X_new, "F"); Zn <- norm(Z_new, "F")
    p_res <- norm(resid, "F")
    d_res <- norm(rho * Dt %*% (Z_new - Z), "F")
    rms_p <- p_res / sqrt(n_edges * n_genes)
    eps_pri  <- sqrt(length(resid)) * tol + tol * max(DX_n, Zn)
    eps_dual <- sqrt(length(X_new)) * tol + tol * norm(rho * Dt %*% u, "F")

    Z <- Z_new; X <- X_new

    if (verbose && iter %% 100L == 0L)
      cat(sprintf("  tree_delta ADMM iter %d: rms_p=%.2e r=%.2e/%.2e xKKT=%.2e rho=%.3g\n",
                  iter, rms_p, p_res, eps_pri, xinfo$kkt_residual, rho))
    if (p_res <= eps_pri && d_res <= eps_dual && xinfo$kkt_residual <= kkt_tol) {
      converged <- TRUE; exit_reason <- "kkt_gate"; break
    }
    if (x_stall_tol > 0 && iter %% chk_every == 0L) {
      if (norm(X - X_ref, "F") / (norm(X_ref, "F") + 1e-12) < x_stall_tol) {
        converged <- TRUE; exit_reason <- "x_stall"; break
      }
      X_ref <- X
    }

    if (adaptive_rho && iter %% 25L == 0L && iter >= 100L) {
      new_rho <- if (p_res > rho_mu * d_res) rho * rho_tau
                 else if (d_res > rho_mu * p_res) rho / rho_tau
                 else rho
      if (new_rho != rho) {
        u <- u * (rho / new_rho); rho <- new_rho
        H <- make_H(rho)
        n_rho_upd <- n_rho_upd + 1L
      }
    }
  }
  if (verbose && !converged)
    cat(sprintf("  tree_delta ADMM: hit max_iter=%d (rms_p=%.2e xKKT=%.2e)\n",
                max_iter, p_res / sqrt(n_edges * n_genes), xinfo$kkt_residual))

  Delta <- tryCatch(solve(T_mat, X), error = function(e) qr.solve(T_mat, X))
  list(X = X, Delta = Delta, T_mat = T_mat,
       residual_norm = norm(Y - Pi %*% X, "F"), iterations = iters,
       rho_final = rho, converged = converged, exit_reason = exit_reason,
       x_kkt_residual = xinfo$kkt_residual, primal_residual = p_res, dual_residual = d_res,
       eps_primal = eps_pri, eps_dual = eps_dual)
}


# ---- GSEA analysis --------------------------------------------------

#' GSEA analysis using fgsea
#'
#' @export
#' @import ggplot2 fgsea ggrepel pheatmap DESeq2
runGSEA <- function(X_optimal,
                    outputDir,
                    treeFile,
                    GSEA_file     = NULL,
                    top_K         = 5,
                    n_permutations = 10000) {

  GSEA_dir <- file.path(outputDir, "GSEA")
  suppressWarnings(dir.create(GSEA_dir))

  X <- as.matrix(X_optimal)
  storage.mode(X) <- "numeric"
  X <- t(X)   # rows = genes, cols = clones

  if (is.null(GSEA_file)) {
    GSEA_file <- system.file("extdata", "h.all.v2024.1.Hs.symbols.gmt.txt",
                              package = "pictographPlus")
  }
  gene_list <- read_GSEA_file(GSEA_file)

  edge_list <- read_tree(treeFile)

  gsea_results_list <- list()

  if (is.null(nrow(edge_list)) || nrow(edge_list) == 1) {
    if (is.null(nrow(edge_list))) {
      sample1 <- as.character(edge_list[1])
      sample2 <- as.character(edge_list[2])
    } else {
      sample1 <- as.character(edge_list[1, 1])
      sample2 <- as.character(edge_list[1, 2])
    }
    gsea_results_list[[1]] <- GSEA_diff(X, sample1, sample2, gene_list,
                                         GSEA_dir, n_permutations, top_K)
    leadingEdgePlot(log2(X + 1), sample1, sample2, GSEA_dir)
  } else {
    edge_list <- apply(edge_list, 2, as.character)
    for (i in seq_len(nrow(edge_list))) {
      sample1 <- edge_list[i, 1]
      sample2 <- edge_list[i, 2]
      gsea_results_list[[i]] <- GSEA_diff(X, sample1, sample2, gene_list,
                                           GSEA_dir, n_permutations, top_K)
      leadingEdgePlot(log2(X + 1), sample1, sample2, GSEA_dir)
    }
  }
}


#' @import purrr
leadingEdgePlot <- function(X, sample1, sample2, GSEA_dir) {
  pdf_filename <- file.path(GSEA_dir,
                             paste0("heatmaps_", sample1, "_", sample2,
                                    "_top30.pdf"))
  pdf(pdf_filename, width = 8, height = 6)

  pathway_data <- read.csv(file.path(GSEA_dir,
                                      paste0("clone", sample2, "_",
                                             sample1, "_GSEA_diff.csv")))

  filtered_pathways <- pathway_data %>%
    filter(padj < 0.05) %>%
    mutate(leadingEdge = strsplit(as.character(leadingEdge), ";")) %>%
    mutate(leadingEdge = map(leadingEdge,
                              ~ .x[seq_len(min(30, length(.x)))])) %>%
    unnest(leadingEdge)

  significant_pathways <- filtered_pathways %>%
    distinct(pathway, NES) %>%
    arrange(desc(NES)) %>%
    pull(pathway)

  for (pathway in significant_pathways) {
    pathway_genes <- filtered_pathways %>%
      filter(pathway == !!pathway) %>%
      pull(leadingEdge) %>%
      unique()

    # Only keep genes present in X
    pathway_genes <- intersect(pathway_genes, rownames(X))
    pathway_X     <- X[pathway_genes, c(sample1, sample2), drop = FALSE]
    if (nrow(pathway_X) == 0) next

    pathway_info <- filtered_pathways %>%
      filter(pathway == !!pathway) %>%
      select(padj, NES) %>%
      distinct()

    pheatmap(
      pathway_X,
      cluster_rows  = FALSE,
      cluster_cols  = TRUE,
      color         = colorRampPalette(c("blue", "white", "red"))(100),
      fontsize_row  = 8,
      fontsize_col  = 8,
      main          = paste0("Pathway: ", pathway, "\n",
                             "padj: ", signif(pathway_info$padj, 3),
                             ", NES: ", signif(pathway_info$NES, 3)),
      border_color  = NA
    )
  }
  dev.off()
}


GSEA_diff <- function(expr_matrix, sample1, sample2, gene_list, GSEA_dir,
                       n_permutations = 10000, n = 5, thresh = 10,
                       minSize = 15, maxSize = 500, eps = 1e-3) {
  # thresh/minSize/maxSize/eps match the documented edge test (response_letter.md
  # R1-M4) and analysis/09_uncertainty_layer's bootstrap scorer: genes expressed
  # (X > thresh) in at least one of the two clones on this edge, tested with the same
  # fgseaMultilevel tolerances, so baseline and bootstrap pathway calls agree.
  keep <- expr_matrix[, sample1] > thresh | expr_matrix[, sample2] > thresh
  expr_matrix <- expr_matrix[keep, , drop = FALSE]

  log2_diff    <- log2((expr_matrix[, sample2] + 1) /
                       (expr_matrix[, sample1] + 1))
  ranked_genes <- sort(log2_diff, decreasing = TRUE)

  gsea_results <- fgseaMultilevel(pathways = gene_list, stats = ranked_genes,
                                   minSize = minSize, maxSize = maxSize, eps = eps)
  gsea_results$Log10padj <- -log10(gsea_results$padj)

  top_up   <- gsea_results %>% arrange(desc(NES)) %>% slice_head(n = n)
  top_down <- gsea_results %>% arrange(NES)        %>% slice_head(n = n)
  gsea_top <- bind_rows(top_up, top_down) %>% filter(padj < 0.05)
  gsea_sig <- gsea_results %>% filter(padj < 0.05)

  sample1N <- if (sample1 == "0") "R" else sample1

  make_bar <- function(df, title) {
    ggplot(df, aes(x = Log10padj, y = reorder(pathway, NES), fill = NES)) +
      geom_bar(stat = "identity") +
      geom_vline(xintercept = -log10(0.05), color = "green",
                 linetype = "dashed", linewidth = 1) +
      scale_fill_gradientn(colors = c("blue", "red"), name = "NES",
                            oob = scales::squish, limits = c(-2, 2)) +
      labs(x = "-log10(p-adj)", y = "Pathway", title = title) +
      theme(axis.text.y  = element_text(size = 24),
            axis.text.x  = element_text(size = 18),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title   = element_text(size = 26))
  }

  title_str <- paste("Clone", sample2, "vs Clone", sample1N)
  ggsave(file.path(GSEA_dir, paste0("clone", sample2, "_", sample1,
                                     "_GSEA_diff.png")),
         plot = make_bar(gsea_top, title_str), width = 13, height = 12)
  ggsave(file.path(GSEA_dir, paste0("clone", sample2, "_", sample1,
                                     "_GSEA_sig.png")),
         plot = make_bar(gsea_sig, title_str), width = 13, height = 12)

  gsea_df             <- as.data.frame(gsea_results)
  gsea_df$leadingEdge <- sapply(gsea_df$leadingEdge, paste, collapse = ";")
  write.csv(gsea_df,
            file.path(GSEA_dir, paste0("clone", sample2, "_", sample1,
                                        "_GSEA_diff.csv")))
  return(gsea_results)
}


read_GSEA_file <- function(GSEA_file) {
  lines     <- readLines(GSEA_file)
  gene_list <- list()
  for (line in lines) {
    elements              <- strsplit(line, "\t")[[1]]
    gene_list[[elements[1]]] <- elements[-c(1, 2)]
  }
  gene_list
}


volcanoplot <- function(data_full, title, n = 5) {
  data_full <- data_full %>%
    mutate(Significance = ifelse(padj < 0.05,
                                 ifelse(NES < 0, "Negative", "Positive"),
                                 "Not Significant"))

  labels <- bind_rows(
    data_full %>% arrange(desc(NES)) %>% slice_head(n = n),
    data_full %>% arrange(NES)       %>% slice_head(n = n)
  )

  ggplot(data_full, aes(x = NES, y = Log10padj, color = Significance)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Negative" = "blue", "Positive" = "red",
                                   "Not Significant" = "gray")) +
    geom_text_repel(data = labels, aes(label = pathway),
                    size = 3.5, max.overlaps = Inf) +
    labs(title = title, x = "Enrichment Score", y = "-log10(P-value)",
         color = "Significant (p < 0.05)") +
    theme_minimal() +
    theme(legend.position = "top")
}
