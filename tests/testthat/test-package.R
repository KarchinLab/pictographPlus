test_that("package loads and exports core entry points", {
  expect_true(requireNamespace("pictographPlus", quietly = TRUE))
  for (fn in c("runPictograph", "runDeconvolution", "runGSEA",
               "runPICTographPlus", "importFiles")) {
    expect_true(exists(fn, where = asNamespace("pictographPlus"),
                       inherits = FALSE),
                info = paste("missing exported function:", fn))
  }
})

test_that("bundled example data files are present", {
  ex_dir <- system.file("extdata", "examples", package = "pictographPlus")
  expect_true(nzchar(ex_dir))
  expect_true(dir.exists(ex_dir))
  expect_true(length(list.files(ex_dir, pattern = "\\.csv$")) > 0)
})

test_that("bundled JAGS model files are present", {
  for (m in c("type1.jags", "type2.jags", "type3.jags",
              "type1_K1.jags", "type2_K1.jags", "type3_K1.jags")) {
    f <- system.file("extdata", m, package = "pictographPlus")
    expect_true(nzchar(f), info = paste("missing JAGS model:", m))
  }
})
