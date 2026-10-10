# Regression tests for the deconvolution solver changes in 1.1.0 (exact X-subproblem, unbounded adaptive rho,
# strict convergence check). See NEWS.md.

kkt_res <- function(H, B, X) { g <- H %*% X - B; max(abs(ifelse(X > 1e-10, g, pmin(g, 0)))) }

test_that(".exact_qp solves the non-negative QP exactly and matches .projected_qp", {
  set.seed(1)
  K <- 6; n <- 200
  A <- matrix(rnorm(K * K), K); H <- crossprod(A) + diag(0.1, K)
  B <- matrix(rnorm(K * n, sd = 3), K, n)
  ex <- pictographPlus:::.exact_qp(H, B, matrix(TRUE, K, n))
  pq <- pictographPlus:::.projected_qp(H, B)
  expect_true(all(ex$X >= 0))
  expect_lt(kkt_res(H, B, ex$X), 1e-8)
  expect_equal(ex$n_fallback, 0)
  expect_lt(max(abs(ex$X - pq$X)), 1e-5)
})

make_problem <- function(K = 5, S = 3, n = 300, seed = 2) {
  set.seed(seed)
  edges <- cbind(0L, seq_len(K - 1L))
  Pi <- matrix(runif(S * K), S, K); Pi <- Pi / rowSums(Pi)
  X0 <- matrix(rexp(K * n, 1 / 50), K, n)
  list(Y = Pi %*% X0 + matrix(rnorm(S * n, sd = 2), S, n), Pi = Pi, edges = edges)
}

test_that("ADMM fits converge at the strict check on an under-determined problem (S < K)", {
  p <- make_problem()
  fits <- list(
    tree_delta  = pictographPlus:::fit_tree_delta_admm(p$Y, p$Pi, p$edges, lambda = 0.01),
    fused_ew    = pictographPlus:::fit_elementwise_fused_lasso_admm(p$Y, p$Pi, p$edges, lambda = 0.01),
    elastic_net = pictographPlus:::fit_elastic_net_tree(p$Y, p$Pi, p$edges, lambda1 = 0.01, lambda2 = 0.01))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_true(f$converged, info = nm)
    expect_identical(f$exit_reason, "kkt_gate", info = nm)
    expect_lte(f$primal_residual, f$eps_primal)
    expect_lte(f$dual_residual, f$eps_dual)
    expect_lte(f$x_kkt_residual, 1e-5)
    expect_true(all(f$X >= 0))
  }
})

test_that("tree_delta fits are deterministic (cold start)", {
  p <- make_problem(seed = 3)
  a <- pictographPlus:::fit_tree_delta_admm(p$Y, p$Pi, p$edges, lambda = 0.01)
  b <- pictographPlus:::fit_tree_delta_admm(p$Y, p$Pi, p$edges, lambda = 0.01)
  expect_identical(a$X, b$X)
  expect_identical(a$iterations, b$iterations)
})
