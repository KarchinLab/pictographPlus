# **PICTographPlus**

## **Overview**


PICTographPlus is a computational tool that integrates bulk DNA and RNA sequencing data to:

1. **Reconstruct Clone-Specific Transcriptomic Profiles**
2. **Infer Tumor Evolution**
3. **Identify Transcriptional Transitions Between Clones**

The tool infers tumor clonal evolution from single or multi-region sequencing data by modeling the uncertainty of mutation cellular fraction (MCF) in small somatic mutations (SSMs) and copy number alterations (CNAs). Using a Bayesian hierarchical model, it assigns SSMs and CNAs to subclones, reconstructing tumor evolutionary trees that adhere to principles of lineage precedence, sum condition, and optional constraints based on sample presence. For deconvolution, PICTographPlus integrates tumor clonal tree structures with clone proportions across samples to resolve bulk gene expression data. It optimizes an objective function that minimizes discrepancies between observed and predicted sample-level gene expression while imposing a smoothness penalty, ensuring that closely related clones display greater gene expression similarity. Lastly, the tool conducts pathway enrichment analysis to identify statistically significant alterations in pathways connecting tumor clones.

### **Core Modules**
- **`runPictograph`** – Tumor evolution inference using genomic data
- **`runDeconvolution`** – Bulk RNA expression deconvolution based on tumor evolution
- **`runGSEA`** – Gene Set Enrichment Analysis (GSEA) for transcriptomic differences between clones

### **Key Features**
- Uses **Bayesian hierarchical modeling** to infer tumor clonal evolution.
- Deconvolves bulk gene expression data using **tumor clonal tree structures** with 7 model variants.
- Performs **pathway enrichment analysis** to highlight significant transcriptomic alterations.

### **Deconvolution Models**

Seven models are available in `runDeconvolution` (and `runPICTographPlus`) via the `model` parameter:

| Model | Default λ | Best for |
|:------|:---------:|:---------|
| `tree_delta` (**default**) | 0.05 | Matched normal sample; best pathway F1 on synthetic edges |
| `plain` | 0.10 | Highest Pearson recovery; fastest solver |
| `adaptive` | 0.50 | General-purpose alternative to plain |
| `adaptive_v2` | 0.50 | Tumor-only mode (no normal reference) |
| `plain_debiased` | 0.50 | Matched normal; best sensitivity |
| `fused_ew` | 0.01 | Highest precision / MCC; requires λ > 0 |
| `elastic_net` | 0.01 | External normal only; requires λ > 0 |

λ values are derived from benchmarking on 320 synthetic pseudo-bulk replicates
(P7 tumour, deep tree, K = 8 clones) using the star-topology as a conservative
topology-agnostic λ selection criterion.

---

## **Bug Fixes (2026-04-16)**

### `deconvolution.R`

1. **Missing row names on `X_optimal`** — The deconvolved expression matrix
   written to `clonal_expression.csv` lacked proper clone-ID row names. When
   read back by `runGSEA`, GSEA comparisons (e.g. clone "0" vs clone "1")
   silently failed because column lookups used 1-based integer names instead
   of clone-ID strings. **Fix:** `rownames(X_optimal) <- colnames(Pi)` is now
   set immediately after solving.

2. **Dead code `tree <- read.csv(treeFile)`** — `runDeconvolution` called
   `read.csv(treeFile)` into an unused variable `tree` before separately
   calling `read_tree(treeFile)`. The CSV read was harmless but wasteful and
   misleading. **Fix:** removed.

3. **Single-edge matrix shape** — `read_tree()` for a single-edge tree
   returned a 2-element vector. The caller did not guarantee matrix shape
   before passing to the model functions. **Fix:** single-edge output is
   coerced to a 1 × 2 integer matrix in `runDeconvolution`.

4. **Gradient-descent solver replaced** — The original `optimize_X`
   gradient-descent loop (100 000 iterations, learning rate 0.01) was slow,
   sensitive to initialisation, and produced approximate solutions. It also
   applied the non-spectral Laplacian from igraph, which differs from the
   normalised form used in benchmarking. **Fix:** replaced with the efficient
   closed-form and ADMM solvers from the validated benchmarking pipeline
   (`shared_scripts/models.R`), now embedded directly in `deconvolution.R`.

5. **`runGSEA` coercion bug** — `mapply(as.matrix(X_optimal), FUN=as.numeric)`
   used `mapply` to coerce matrix elements, which is non-standard and fails on
   certain matrix classes. **Fix:** replaced with `as.matrix()` +
   `storage.mode(X) <- "numeric"`.

6. **`size` deprecated in `geom_vline`** — `ggplot2` ≥ 3.4 deprecates `size`
   for line geoms in favour of `linewidth`. **Fix:** updated to `linewidth = 1`.

7. **`leadingEdgePlot` hard-coded gene filter** — `1:min(30, length(.x))`
   used base R indexing without `seq_len`, which can produce unexpected results
   for zero-length vectors. **Fix:** replaced with `seq_len(min(30, length(.x)))`.

---

## **Installation**

### **Install JAGS (Required for Bayesian Analysis)**
JAGS must be installed separately. Download it from: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)

### **Install PICTographPlus**
Run the following command in R:

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("KarchinLab/pictographPlus", build_vignettes = TRUE)
```

### **Package versions during development**
PICTographPlus was developed under R (4.4.2). All package versions during development can be found at [installed_packages.csv](https://github.com/KarchinLab/pictographPlus/blob/working/installed_packages.csv)

--- 

## **Access Tutorial**

Detailed tutorial can be accessed through [vignette](https://github.com/KarchinLab/pictographPlus/tree/working/vignettes).
```r
library(pictographPlus)
vignette("pictographPlus", package = "pictographPlus")
```

