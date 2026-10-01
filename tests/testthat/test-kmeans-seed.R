# Regression tests for the k-means step in check_sample_CNA (allele-specific copy number from
# heterozygous SNVs). Before 2026-09-30 it called kmeans() unseeded with one start, so the processed
# copy-number input (and the tree) could change between identical runs. See preprocessing.R.

# A segment whose tumor VAFs are clearly bimodal or trimodal, so check_sample_CNA reaches kmeans.
make_inputs <- function(dir, modes = c(0.2, 0.8), n_snp = 60, seed = 1) {
  set.seed(seed)
  depth <- 200
  vaf <- sample(modes, n_snp, replace = TRUE) + rnorm(n_snp, 0, 0.03)
  alt <- round(depth * pmin(pmax(vaf, 0.01), 0.99))
  snp <- data.frame(chrom = "chr1", position = seq(1000, by = 1000, length.out = n_snp),
                    germline_ref = 50, germline_alt = 50,
                    sample1_ref = depth - alt, sample1_alt = alt)
  snp_file <- file.path(dir, "snp.csv")
  utils::write.csv(snp, snp_file, row.names = FALSE)
  cna <- tibble::tibble(sample = "sample1", chrom = "chr1", start = 1, end = 1e6, tcn = 3.2)
  list(cna = cna, snp_file = snp_file)
}

run_check <- function(inp, dir, ...) {
  suppressMessages(pictographPlus:::check_sample_CNA(inp$cna, dir, inp$snp_file, LOH = FALSE, ...))
}

test_that("kmeans_seeded is deterministic and leaves the global RNG untouched", {
  x <- c(rnorm(30, 0.2, 0.05), rnorm(30, 0.5, 0.05), rnorm(30, 0.8, 0.05))
  set.seed(42); before <- .Random.seed
  a <- pictographPlus:::kmeans_seeded(x, centers = 3, seed = 7)
  expect_identical(.Random.seed, before)
  set.seed(999)
  b <- pictographPlus:::kmeans_seeded(x, centers = 3, seed = 7)
  expect_identical(a$cluster, b$cluster)
  expect_identical(a$centers, b$centers)
})

test_that("kmeans_seeded does not depend on the session's RNGkind", {
  x <- c(rnorm(30, 0.2, 0.05), rnorm(30, 0.8, 0.05))
  a <- pictographPlus:::kmeans_seeded(x, centers = 2, seed = 7)
  old <- RNGkind("L'Ecuyer-CMRG"); on.exit(RNGkind(old[1], old[2], old[3]))
  b <- pictographPlus:::kmeans_seeded(x, centers = 2, seed = 7)
  expect_identical(a$cluster, b$cluster)
})

test_that("check_sample_CNA output is identical across global RNG states", {
  for (modes in list(c(0.2, 0.8), c(0.15, 0.5, 0.85))) {
    dir <- withr::local_tempdir()
    inp <- make_inputs(dir, modes)
    # 30 global seeds: with the old unseeded single-start kmeans, the trimodal input gave 3
    # different tcn_alt values over these seeds, so this would have failed before the fix.
    set.seed(123); state <- .Random.seed
    r1 <- run_check(inp, dir)
    expect_identical(.Random.seed, state)         # global stream not consumed
    expect_true(r1$tcn_alt[1] > 0)                # the kmeans branch was actually used
    for (s in 1:30) { set.seed(s); expect_identical(run_check(inp, dir), r1) }
  }
})

test_that("the same preprocessing_seed reproduces; nstart = 50 makes the result start-independent", {
  dir <- withr::local_tempdir()
  inp <- make_inputs(dir, c(0.15, 0.5, 0.85))
  ref <- run_check(inp, dir, preprocessing_seed = 123)
  for (s in c(1, 2, 3, 456, 2027)) {
    expect_identical(run_check(inp, dir, preprocessing_seed = s)[, c("tcn_alt", "tcn_ref", "to_keep")],
                     ref[, c("tcn_alt", "tcn_ref", "to_keep")])
  }
})
