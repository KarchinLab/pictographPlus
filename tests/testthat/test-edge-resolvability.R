# edgeResolvability() (development version): an edge between clones with their own pure samples is fully resolvable
# (noise multiplier 1); clones whose proportions move together are not.
test_that("edgeResolvability separates resolvable from confounded edges", {
  d <- withr::local_tempdir()
  write.csv(data.frame(edge = c("root->1", "1->2", "1->3"), parent = c("root", "1", "1"), child = c("1", "2", "3")),
            file.path(d, "tree.csv"), row.names = FALSE, quote = FALSE)
  pr <- data.frame(s1 = c(1, 0, 0), s2 = c(0, 1, 0), s3 = c(0, 0, 1), row.names = 1:3)
  write.csv(pr, file.path(d, "prop.csv"))
  r <- edgeResolvability(file.path(d, "prop.csv"), file.path(d, "tree.csv"), normal_samples = "n1")
  expect_equal(r$noise_multiplier, c(1, 1, 1)); expect_equal(r$estimable_fraction, c(1, 1, 1)); expect_true(all(r$resolvable))
  pr2 <- data.frame(s1 = c(0.5, 0.5, 0), s2 = c(0.25, 0.25, 0.5), row.names = 1:3)          # clones 1 and 2 always equal
  write.csv(pr2, file.path(d, "prop2.csv"))
  r2 <- edgeResolvability(file.path(d, "prop2.csv"), file.path(d, "tree.csv"))
  expect_lt(r2$estimable_fraction[r2$edge == "clone2_1"], 0.5); expect_false(r2$resolvable[r2$edge == "clone2_1"])
})
