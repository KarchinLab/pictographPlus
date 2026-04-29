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
