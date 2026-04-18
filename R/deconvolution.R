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
#' @param lambda Regularisation strength; default 0.01 (elastic_net star-best).
#' @param use_star_tree If TRUE (default), ignore \code{treeFile} and use a
#'   star topology (root directly connected to all clones). Recommended when
#'   the inferred tree topology does not improve performance. When FALSE,
#'   \code{treeFile} must be supplied.
#' @param model Deconvolution model. One of:
#'   \itemize{
#'     \item \code{"elastic_net"} (default) — L2 Laplacian + L1 fused LASSO
#'       via ADMM. Best synthetic-edge F1 (0.347) and sensitivity (0.368) in
#'       with-normal mode; requires lambda > 0.
#'     \item \code{"tree_delta"} — tree-delta parameterisation with group-L2
#'       ADMM. Near-tied F1 (0.339) in with-normal; best F1 in with-extnorm
#'       mode. Use when interpretability of tree-structured penalty is valued.
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
                             model          = "elastic_net",
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
    purityData   <- read.csv(purityFile)
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
  proportionDF <- proportionDF[order(rownames(proportionDF)), , drop = FALSE]

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


# ---- Model 5: Element-wise fused LASSO via ADMM --------------------

fit_elementwise_fused_lasso_admm <- function(Y, Pi, edges, lambda = 0.01,
                                              rho = 1.0, max_iter = 5000,
                                              tol = 1e-3, ridge = 1e-8,
                                              verbose = FALSE) {
  K       <- max(edges) + 1L
  n_genes <- ncol(Y)
  n_edges <- nrow(edges)
  D       <- build_incidence(edges, K)

  A     <- 2.0 * t(Pi) %*% Pi + rho * t(D) %*% D + ridge * diag(K)
  A_inv <- tryCatch(solve(A),
                    error = function(e) qr.solve(A, diag(K)))

  X <- matrix(0.0, K, n_genes)
  Z <- matrix(0.0, n_edges, n_genes)
  u <- matrix(0.0, n_edges, n_genes)
  tau   <- lambda / rho
  iters <- 0L

  for (iter in seq_len(max_iter)) {
    iters  <- iter
    X_new  <- pmax(A_inv %*% (2.0 * t(Pi) %*% Y + rho * t(D) %*% (Z - u)), 0)
    V      <- D %*% X_new + u
    Z_new  <- sign(V) * pmax(abs(V) - tau, 0)
    resid  <- D %*% X_new - Z_new
    u      <- u + resid

    DX_n   <- norm(D %*% X_new, "F"); Zn <- norm(Z_new, "F")
    d_res  <- norm(rho * t(D) %*% (Z_new - Z), "F")
    d_scl  <- norm(rho * t(D) %*% u, "F") + 1e-10
    rel_p  <- norm(resid, "F") / (max(DX_n, Zn) + 1e-10)
    rel_d  <- d_res / d_scl
    Z <- Z_new; X <- X_new

    if (verbose && iter %% 100 == 0)
      cat(sprintf("  fused_ew ADMM iter %d: rel_p=%.2e rel_d=%.2e\n",
                  iter, rel_p, rel_d))
    if (rel_p < tol && rel_d < tol) break
  }
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"), iterations = iters)
}


# ---- Model 6: Elastic net (L2 + L1) via ADMM -----------------------

fit_elastic_net_tree <- function(Y, Pi, edges, lambda1 = 0.01, lambda2 = 0.01,
                                  ridge = 1e-8, normalize = "spectral",
                                  max_iter = 5000, tol = 1e-3,
                                  verbose = FALSE) {
  K       <- max(edges) + 1L
  n_genes <- ncol(Y)
  n_edges <- nrow(edges)
  L       <- build_laplacian(edges, K, normalize = normalize)
  D       <- build_incidence(edges, K)

  rho   <- 1.0
  A     <- 2.0 * t(Pi) %*% Pi + lambda1 * L + rho * t(D) %*% D + ridge * diag(K)
  A_inv <- tryCatch(solve(A),
                    error = function(e) qr.solve(A, diag(K)))

  X <- matrix(0.0, K, n_genes)
  Z <- matrix(0.0, n_edges, n_genes)
  u <- matrix(0.0, n_edges, n_genes)
  tau   <- lambda2 / rho
  iters <- 0L

  for (iter in seq_len(max_iter)) {
    iters  <- iter
    X_new  <- pmax(A_inv %*% (2.0 * t(Pi) %*% Y + rho * t(D) %*% (Z - u)), 0)
    V      <- D %*% X_new + u
    Z_new  <- sign(V) * pmax(abs(V) - tau, 0)
    resid  <- D %*% X_new - Z_new
    u      <- u + resid

    DX_n  <- norm(D %*% X_new, "F"); Zn <- norm(Z_new, "F")
    d_res <- norm(rho * t(D) %*% (Z_new - Z), "F")
    d_scl <- norm(rho * t(D) %*% u, "F") + 1e-10
    rel_p <- norm(resid, "F") / (max(DX_n, Zn) + 1e-10)
    rel_d <- d_res / d_scl
    Z <- Z_new; X <- X_new

    if (verbose && iter %% 100 == 0)
      cat(sprintf("  elastic_net ADMM iter %d: rel_p=%.2e rel_d=%.2e\n",
                  iter, rel_p, rel_d))
    if (rel_p < tol && rel_d < tol) break
  }
  list(X = X, residual_norm = norm(Y - Pi %*% X, "F"), iterations = iters)
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
                                 rho = 1.0, max_iter = 3000, tol = 1e-4,
                                 ridge = 1e-8, verbose = FALSE) {
  K       <- max(edges) + 1L
  n_genes <- ncol(Y)
  n_edges <- nrow(edges)

  T_mat    <- build_path_matrix(edges, K)
  Pi_tilde <- Pi %*% T_mat
  PtP      <- t(Pi_tilde) %*% Pi_tilde
  PtY      <- t(Pi_tilde) %*% Y

  if (lambda == 0) {
    A     <- PtP + ridge * diag(K)
    Delta <- tryCatch(solve(A, PtY), error = function(e) qr.solve(A, PtY))
    X     <- pmax(T_mat %*% Delta, 0)
    return(list(X = X, Delta = Delta, T_mat = T_mat,
                residual_norm = norm(Y - Pi %*% X, "F"), iterations = 0L))
  }

  C    <- diag(c(0.0, rep(1.0, n_edges)))
  A    <- 2.0 * PtP + rho * C + ridge * diag(K)
  Ainv <- tryCatch(solve(A),
                   error = function(e) qr.solve(A, diag(K)))

  Delta <- matrix(0.0, K, n_genes)
  V     <- matrix(0.0, n_edges, n_genes)
  U     <- matrix(0.0, n_edges, n_genes)
  tau   <- lambda / rho
  iters <- 0L

  for (iter in seq_len(max_iter)) {
    iters <- iter
    RHS              <- 2.0 * PtY
    RHS[2L:K, ]     <- RHS[2L:K, ] + rho * (V - U)
    Delta_new        <- Ainv %*% RHS

    delta_e <- Delta_new[2L:K, , drop = FALSE]
    W       <- delta_e + U
    V_new   <- matrix(0.0, n_edges, n_genes)
    for (e in seq_len(n_edges)) {
      wn <- sqrt(sum(W[e, ]^2))
      if (wn > tau) V_new[e, ] <- W[e, ] * (1.0 - tau / wn)
    }

    resid <- delta_e - V_new
    U_new <- U + resid
    rms_p <- norm(resid, "F") / sqrt(n_edges * n_genes)
    d_res <- norm(rho * (V_new - V), "F")

    V <- V_new; U <- U_new; Delta <- Delta_new

    if (verbose && iter %% 100L == 0L)
      cat(sprintf("  tree_delta ADMM iter %d: rms_p=%.2e dual=%.2e\n",
                  iter, rms_p, d_res))
    if (rms_p < tol) break
  }

  X <- pmax(T_mat %*% Delta, 0)
  list(X = X, Delta = Delta, T_mat = T_mat,
       residual_norm = norm(Y - Pi %*% X, "F"), iterations = iters)
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
                       n_permutations = 10000, n = 5) {
  log2_diff    <- log2((expr_matrix[, sample2] + 1) /
                       (expr_matrix[, sample1] + 1))
  ranked_genes <- sort(log2_diff, decreasing = TRUE)

  gsea_results <- fgseaMultilevel(pathways = gene_list, stats = ranked_genes)
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
