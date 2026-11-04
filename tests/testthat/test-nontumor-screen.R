# screenNonTumor() (1.2.0): flags a tumor clone that carries non-tumor marker genes, and pathway calls whose leading
# edge is made of non-tumor-specific genes.
test_that("screenNonTumor flags a clone carrying hepatocyte genes and a stroma-driven pathway call", {
  ref <- read.csv(system.file("extdata", "nontumor_genes.csv", package = "pictographPlus"))
  liver <- c("ALB", "APOA1", "FGB", "SERPINA1", "TF"); strom <- head(ref$gene[ref$group == "stroma"], 10)
  genes <- unique(c(liver, strom, paste0("TUMOR", 1:40)))
  X <- matrix(100, 3, length(genes), dimnames = list(c("0", "1", "2"), genes))
  X["1", liver] <- 0; X["2", liver] <- 5000; X["0", liver] <- 200
  d <- withr::local_tempdir()
  write.csv(data.frame(pathway = c("STROMA_SET", "TUMOR_SET", "NS_SET"), padj = c(0.01, 0.01, 0.5), NES = c(2, -2, 1),
                       leadingEdge = c(paste(strom, collapse = ";"), paste0("TUMOR", 1:10, collapse = ";"), "TUMOR1")),
            file.path(d, "clone2_1_GSEA_diff.csv"), row.names = FALSE)
  out <- suppressMessages(capture.output(r <- screenNonTumor(X, GSEA_dir = d, sites = "liver", outputDir = d)))
  hep <- r$markers[r$markers$cell_type == "hepatocyte", ]
  expect_equal(hep$flag, c(FALSE, FALSE, TRUE))                       # clones 0, 1, 2
  expect_equal(nrow(r$pathways), 2L)                                  # the non-significant call is dropped
  expect_equal(r$pathways$flag, c(TRUE, FALSE))
  expect_equal(r$pathways$frac_stroma[1], 1)
  expect_true(file.exists(file.path(d, "nontumor_screen_pathways.csv")))
  expect_error(screenNonTumor(X, sites = "spleen"), "unknown site")
  expect_error(screenNonTumor(X[-1, ]), "row \"0\"")
})
