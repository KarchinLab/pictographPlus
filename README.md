# **PICTographPlus**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19896731.svg)](https://doi.org/10.5281/zenodo.19896731)

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

#### Model Selection Quick-Reference

| Scenario | Recommended model | Evidence (companion study) |
|:---------|:------------------|:---------------------------|
| Default: matched normal sample available, tree uncertain | `tree_delta` on the star graph (λ = 0.01; package default) | Higher clone-specific recovery than same-input NNLS in 7 of 8 wellDR-seq patients. Aggregate recovery is similar to minimum-norm least squares, so the gain over NNLS is not attributed to the graph penalty alone. |
| Tree confidently known (e.g., stable across MCMC chains) | Consider `tree_delta` or `fused_ew` on the inferred DNA tree (`use_star_tree = FALSE`) | Higher recovery with the benchmark reference tree in some comparisons. |
| As many samples as clones or more, noisy data | `elastic_net` or `plain` (star graph) | More shrinkage stabilized the estimate under count noise when there were at least as many samples as clones. |
| Clone compositions similar across samples (ill-conditioned proportions) | `elastic_net`; run `edgeResolvability()` first | Most robust in the operating-limits simulation when the proportion matrix was ill-conditioned. |
| No matched normal sample | Any model; do not interpret normal-to-tumor edges | Normal-to-tumor edge F1 was near zero without a normal anchor. |
| Only an unmatched (population-average) normal available | Assess reference compatibility before interpreting normal-to-tumor shifts | The tested atlas-derived normal reduced recovery at low purity; this does not establish the performance of other external references. |
| Any tumor-to-tumor claim | Report by `edgeResolvability()` tier; screen with `screenNonTumor()` | In organoid calibration, signed precision was stratified by tier and recall was low; tier precision is benchmark-specific. |

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
PICTographPlus was developed under R (4.4.2). All package versions during development can be found at [installed_packages.csv](https://github.com/KarchinLab/pictographPlus/blob/main/installed_packages.csv)

---

## **Access Tutorial**

Detailed tutorial can be accessed through [vignette](https://github.com/KarchinLab/pictographPlus/tree/main/vignettes).
```r
library(pictographPlus)
vignette("pictographPlus", package = "pictographPlus")
```

---

## **Citation**

If you use PICTographPlus in your research, please cite the archived software release (the DOI below always resolves to the latest archived version):

> Lai J, Yang Y, Karchin R (2026). *pictographPlus: Reconstructing Clone-Resolved Transcriptional Programs from Bulk Tumor Sequencing.* R package version 1.2.0. Zenodo. https://doi.org/10.5281/zenodo.19896731

You can also retrieve a `bibentry` from R with:

```r
citation("pictographPlus")
```

## **Data**

Processed data generated in the companion study — clone-level expression matrices, inferred clonal trees, and edge-level GSEA results — are available on Mendeley Data at https://doi.org/10.17632/cv66sgfcfn.1.

The benchmarks and applications in the manuscript use the following public datasets:

| Dataset | Accession |
|---|---|
| wellDR-seq scDNA/scRNA co-profile (benchmarking) | GEO: GSE261713 |
| IPMN WES + RNA-seq | dbGaP: phs002225.v3.p1 |
| TRACERx NSCLC WES | EGA: EGAD00001009825 |
| TRACERx NSCLC RNA-seq | EGA: EGAD00001009862 |
| PanCuRx PDAC WGS | EGA: EGAD00001004551 |
| PanCuRx PDAC RNA-seq | EGA: EGAD00001004548 |

All datasets were used in accordance with the respective data access agreements.

