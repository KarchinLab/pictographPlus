# runGSEA()/GSEA_diff() seed (1.1.0): identical pathway tables across runs, and the caller's RNG state untouched.

make_expr <- function(seed = 4, n = 2000) {
  set.seed(seed)
  genes <- sprintf("G%04d", seq_len(n))
  m <- cbind(`0` = rexp(n, 1 / 100), `1` = rexp(n, 1 / 100))
  rownames(m) <- genes
  sets <- lapply(1:20, function(i) sample(genes, 40)); names(sets) <- sprintf("SET_%02d", 1:20)
  m[sets$SET_01, "1"] <- m[sets$SET_01, "1"] * 4
  list(m = m, sets = sets)
}

test_that("seeded GSEA_diff is reproducible and leaves the global RNG untouched", {
  skip_if_not_installed("fgsea")
  e <- make_expr()
  d1 <- withr::local_tempdir(); d2 <- withr::local_tempdir()
  set.seed(42); before <- .Random.seed
  r1 <- suppressWarnings(pictographPlus:::GSEA_diff(e$m, "0", "1", e$sets, d1, seed = 1))
  expect_identical(.Random.seed, before)
  set.seed(7)
  r2 <- suppressWarnings(pictographPlus:::GSEA_diff(e$m, "0", "1", e$sets, d2, seed = 1))
  expect_identical(as.data.frame(r1)[, c("pathway", "pval", "padj", "NES")],
                   as.data.frame(r2)[, c("pathway", "pval", "padj", "NES")])
})
