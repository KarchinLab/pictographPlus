# ---- Non-tumor contamination screen ----------------------------------
#
# The deconvolution represents all non-tumor cells of a patient with ONE expression profile (clone "0"), because DNA
# fixes how much of each sample is non-tumor (1 - purity) but not which cell types those are. When the make-up of
# the non-tumor part differs between samples in step with a clone (e.g. a clone that dominates a metastasis in
# another organ, or a region rich in immune cells or stroma), the difference can be assigned to that clone. This
# screen looks for the fingerprint of that error, using public references only (no matched normal RNA per site).

#' Screen clone expression and edge pathway calls for non-tumor contamination
#'
#' Two checks on the output of \code{runDeconvolution()} and \code{runGSEA()}:
#' \enumerate{
#'   \item \strong{Markers, per clone.} Mean expression of canonical marker genes of immune, stromal, vascular and
#'     (for the \code{sites} given) site-parenchymal cell types in every clone. A tumor clone carrying such genes at
#'     a level comparable to the non-tumor clone "0" is absorbing non-tumor cells.
#'   \item \strong{Pathways, per edge.} For every significant pathway on every edge, the fraction of its
#'     leading-edge genes that are specific to non-tumor cell types (Human Protein Atlas single-cell types: >= 10
#'     nCPM in the cell type and >= 5x every epithelial/tumor-like cell type; cell-cycle genes never listed; see
#'     \code{data-raw/nontumor_reference.R}). Calls above \code{flag_fraction} are flagged: they may reflect
#'     differences in non-tumor content between samples rather than clone biology.
#' }
#' The screen flags; it does not correct. Removing these genes before the pathway test does not remove the
#' contamination, because non-tumor cells change the expression of many genes beyond their markers.
#'
#' @param X_optimal Clone x gene matrix from \code{runDeconvolution()} (row "0" = non-tumor), or the path of its
#'   \code{clonal_expression.csv}.
#' @param GSEA_dir Directory with \code{runGSEA()} outputs (\code{clone<child>_<parent>_GSEA_diff.csv}); NULL skips
#'   the pathway check.
#' @param sites Sites sampled besides the primary organ's tumor, adding their parenchymal genes: any of
#'   "liver", "lung", "pancreas", "brain", "kidney", "adrenal", "bone_marrow". Lymph nodes are covered by the
#'   immune genes, which are always included with stroma and vascular genes.
#' @param outputDir If given, writes \code{nontumor_screen_markers.csv} and \code{nontumor_screen_pathways.csv}.
#' @param padj_cutoff Significance cutoff for pathway calls (default 0.05).
#' @param flag_fraction Leading-edge fraction of non-tumor genes above which a call is flagged (default 0.5).
#' @param marker_min,marker_ratio A tumor clone's marker set is flagged when its mean expression is at least
#'   \code{marker_min} (normalized counts) and at least \code{marker_ratio} times the level in clone "0" (by
#'   default: a tumor clone carrying more of a non-tumor cell type's markers than the non-tumor clone itself).
#' @return A list with \code{markers} (cell type x clone, with flags) and \code{pathways} (one row per
#'   significant call, with the non-tumor fraction overall and per group, and a flag), invisibly; a summary is
#'   printed.
#' @export
screenNonTumor <- function(X_optimal, GSEA_dir = NULL, sites = character(0), outputDir = NULL,
                           padj_cutoff = 0.05, flag_fraction = 0.5, marker_min = 100, marker_ratio = 1) {
  if (is.character(X_optimal)) X_optimal <- read.csv(X_optimal, row.names = 1, check.names = FALSE)
  X <- as.matrix(X_optimal); storage.mode(X) <- "numeric"
  if (!"0" %in% rownames(X)) stop("X_optimal must have a row \"0\" (the non-tumor clone), as returned by runDeconvolution()")
  ok_sites <- c("liver", "lung", "pancreas", "brain", "kidney", "adrenal", "bone_marrow")
  bad <- setdiff(sites, ok_sites)
  if (length(bad)) stop("unknown site(s): ", paste(bad, collapse = ", "), "; use any of ", paste(ok_sites, collapse = ", "))
  groups <- c("immune", "stroma", "vascular", sites)

  # 1. markers per clone
  mk <- read.csv(system.file("extdata", "nontumor_markers.csv", package = "pictographPlus"), stringsAsFactors = FALSE)
  mk <- mk[mk$group %in% groups & mk$gene %in% colnames(X), ]
  markers <- do.call(rbind, lapply(split(mk, factor(mk$cell_type, levels = unique(mk$cell_type))), function(m) {
    lev <- rowMeans(X[, m$gene, drop = FALSE])
    data.frame(cell_type = m$cell_type[1], group = m$group[1], genes = paste(m$gene, collapse = ";"),
               clone = rownames(X), mean_expression = round(lev, 1),
               ratio_to_nontumor = if (lev[["0"]] >= 1) round(lev / lev[["0"]], 3) else NA_real_,
               flag = rownames(X) != "0" & lev >= marker_min & lev >= marker_ratio * lev[["0"]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))

  # 2. pathway calls per edge
  pathways <- NULL
  if (!is.null(GSEA_dir)) {
    ref <- read.csv(system.file("extdata", "nontumor_genes.csv", package = "pictographPlus"), stringsAsFactors = FALSE)
    ref <- ref[ref$group %in% groups, ]; sets <- split(ref$gene, ref$group); all_nt <- unique(ref$gene)
    files <- list.files(GSEA_dir, pattern = "^clone.+_GSEA_diff\\.csv$", full.names = TRUE)
    if (!length(files)) stop("no clone*_GSEA_diff.csv files in ", GSEA_dir)
    pathways <- do.call(rbind, lapply(files, function(f) {
      g <- read.csv(f, stringsAsFactors = FALSE); g <- g[!is.na(g$padj) & g$padj < padj_cutoff, , drop = FALSE]
      if (!nrow(g)) return(NULL)
      le <- strsplit(as.character(g$leadingEdge), "[;,[:space:]]+")
      fr <- sapply(le, function(v) mean(v %in% all_nt))
      byg <- sapply(names(sets), function(s) sapply(le, function(v) round(mean(v %in% sets[[s]]), 3)))
      byg <- matrix(byg, nrow = nrow(g), dimnames = list(NULL, paste0("frac_", names(sets))))
      data.frame(edge = sub("_GSEA_diff\\.csv$", "", basename(f)), pathway = g$pathway,
                 direction = ifelse(g$NES > 0, "UP", "DOWN"), NES = round(g$NES, 3), padj = signif(g$padj, 3),
                 n_leading = lengths(le), frac_nontumor = round(fr, 3), byg, flag = fr > flag_fraction,
                 stringsAsFactors = FALSE)
    }))
  }

  if (!is.null(outputDir)) {
    write.csv(markers, file.path(outputDir, "nontumor_screen_markers.csv"), row.names = FALSE)
    if (!is.null(pathways)) write.csv(pathways, file.path(outputDir, "nontumor_screen_pathways.csv"), row.names = FALSE)
  }
  fm <- markers[markers$flag, ]
  cat(sprintf("Non-tumor screen (groups: %s)\n", paste(groups, collapse = ", ")))
  if (nrow(fm)) {
    cat("Tumor clones carrying more non-tumor marker expression than the non-tumor clone \"0\":\n")
    for (k in unique(fm$clone)) cat(sprintf("  clone %s: %s\n", k, paste(fm$cell_type[fm$clone == k], collapse = ", ")))
  } else cat("No tumor clone carries non-tumor marker genes above the thresholds.\n")
  if (!is.null(pathways)) cat(sprintf("Pathway calls flagged (> %.0f%% non-tumor leading-edge genes): %d of %d\n",
                                      100 * flag_fraction, sum(pathways$flag), nrow(pathways)))
  invisible(list(markers = markers, pathways = pathways))
}
