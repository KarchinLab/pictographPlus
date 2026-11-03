# pictographPlus 1.2.0

## Non-tumor contamination screen (`screenNonTumor()`)

* New function `screenNonTumor()`. The deconvolution represents all non-tumor cells with one profile (clone "0"),
  so non-tumor content that differs between samples in step with a clone can be assigned to that clone. The screen
  reports, per tumor clone, marker genes of immune, stromal, vascular and (optional `sites`) site-parenchymal cell
  types, and, per edge pathway call, the fraction of leading-edge genes specific to non-tumor cell types (Human
  Protein Atlas single-cell types; `inst/extdata/nontumor_genes.csv`, built by `data-raw/nontumor_reference.R`).
  It flags; it does not correct.

## Edge resolvability (`edgeResolvability()`)

* New function `edgeResolvability()`: from the proportions alone (no expression data), how well the samples
  determine each tree edge's clone contrast. It reports the estimable fraction of the contrast (below 1 when part of
  the edge's difference is set by regularization, not data) and a noise multiplier (1 when each clone has a pure
  sample; larger as the two clones' proportions move together).
* Calibrated on two benchmarks with ground truth (wellDR breast cancer; Roerink colorectal organoids): signed
  precision of tumor-to-tumor pathway calls rises with the estimable fraction (Roerink 0.41 / 0.67 / 0.76 for
  < 0.8 / 0.8-0.95 / >= 0.95; wellDR 0.06 / 0.18 / 0.49), which sets the default `min_estimable = 0.95`. The help
  page gives both tables.
* The default noise cut-off is now `max_noise = 100` (was `Inf`): on the Roerink organoids, among edges >= 0.95
  estimable, precision fell from 0.89 / 0.78 to 0.51 above a multiplier of 100. The flag now matches the tiers used
  in the paper. Use `max_noise = Inf` for the previous behavior.
* `screenNonTumor()` help page: benchmark precision of flagged calls, and why stroma needs the per-clone marker
  check rather than the pathway flag.

## Deconvolution solver (`runDeconvolution()`, models `tree_delta`, `fused_ew`, `elastic_net`)

* The ADMM fits now stop only when the full convergence check passes (primal and dual residuals within
  `tol = 1e-6`, and the X-subproblem's KKT residual within `kkt_tol = 1e-5`), with up to `max_iter = 200000`
  iterations (most fits need far fewer: every benchmark fit converged within 60,000; a few under-determined
  real-case fits need just over 60,000). Previously `tol = 1e-4` and `max_iter = 15000`, which on problems with fewer samples than clones
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
* The edge test itself is unchanged: genes expressed (> 10) in ANY clone, fgseaMultilevel's default eps (1e-50),
  as in the benchmark behind the paper's accuracy numbers. (An interim change on the development branch -- genes
  > 10 in the edge's two clones only, eps 1e-3 -- lowered benchmark tumor-to-tumor edge F1 from 0.25 to 0.05 and
  was reverted.)

## Preprocessing

* The allele-specific copy-number k-means step is seeded (`kmeans_seeded()`), so identical inputs give identical
  copy-number calls and trees.

# pictographPlus 1.1.1

* Submission-prep release for Nature Communications.
* Add `LICENSE` file (MIT) and proper `Authors@R` block.
* Restore missing `rjags` and `epiR` Imports; drop `tidyverse` meta-package import.
* Move `getPileUp.py` to `inst/scripts/` and update vignette accordingly.
* Add `tests/testthat` smoke tests, GitHub Actions `R-CMD-check` workflow, `NEWS.md`, and `inst/CITATION`.

# pictographPlus 1.1.0

* Added expanded set of deconvolution model variants (`elastic_net`, `tree_delta`, `adaptive`, `adaptive_v2`, `plain`) with model-selection guidance in the README.
* Improved single-clone handling in `plotSubclone` plotting routines.
* Updates to `runDeconvolution`, `MCMC-main`, `MCMC-clustering`, and `MCMC_plot_clusters` for stability and edge-case handling.

# pictographPlus 1.0.0

* First public release accompanying the manuscript submission.
* Three core modules: `runPictograph` (Bayesian tumor evolution inference), `runDeconvolution` (clone-resolved bulk RNA deconvolution), `runGSEA` (pathway enrichment between clones).
* Vignette with end-to-end example.

# pictographPlus 0.2.0

* Pre-release internal version.

# pictographPlus 0.1.0

* Initial development version.
