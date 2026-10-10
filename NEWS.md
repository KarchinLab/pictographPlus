# pictographPlus 1.1.0

## Deconvolution solver (`runDeconvolution()`, models `tree_delta`, `fused_ew`, `elastic_net`)

* The ADMM fits now stop only when the full convergence check passes (primal and dual residuals within
  `tol = 1e-6`, and the X-subproblem's KKT residual within `kkt_tol = 1e-5`), with up to `max_iter = 60000`
  iterations. Previously `tol = 1e-4` and `max_iter = 15000`, which on problems with fewer samples than clones
  could stop before the optimum. The X-plateau early exit is now an argument, `x_stall_tol`, off by default.
* The X-subproblem (a non-negative least-squares problem per gene) is solved exactly by block principal pivoting
  (Kim & Park 2011), warm-started between iterations, instead of projected gradient descent. Results are the same;
  each iteration is 10-25x faster.
* Adaptive rho (residual balancing) is no longer bounded to `[rho0/64, rho0*64]`, frozen after half of `max_iter`,
  or limited to 40 changes. Those limits kept under-determined `tree_delta` fits from converging.
* `runDeconvolution()` writes `deconvolution_fit_info.csv` (converged, exit reason, iterations, KKT residual) next
  to `clonal_expression.csv`, and warns when a fit reaches `max_iter`.
* `tree_delta` is the default model.

## Pathway analysis (`runGSEA()`)

* New argument `seed` (default 1): each edge's fgsea test is seeded, so pathway calls are reproducible (unseeded,
  about 3-6% of padj < 0.05 calls changed between reruns on identical input). The caller's random-number state is
  left unchanged. `seed = NULL` restores the old behaviour.
* The edge test keeps genes expressed (> 10) in either clone of that edge, with fgseaMultilevel minSize 15,
  maxSize 500, eps 1e-3.

## Preprocessing

* The allele-specific copy-number k-means step is seeded (`kmeans_seeded()`), so identical inputs give identical
  copy-number calls and trees.

# pictographPlus 1.0.0

* Version used for the original submission.
