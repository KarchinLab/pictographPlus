# PICTographPlus
PICTographPlus is a computational tool that integrates bulk DNA and RNA sequencing to reconstruct clone-specific transcriptomic profiles and uncover transcriptional transitions between tumor clones, providing a robust framework for studying cancer evolution. It comprises three core models: 1) tumor evolution inference using genomic data (runPictograph), 2) bulk RNA expression deconvolution based on tumor evolution (runDeconvolution), and 3) gene set enrichment analysis (GSEA) to examine transcriptomic differences between clones (runGSEA).

The tool infers tumor clonal evolution from single or multi-region sequencing data by modeling the uncertainty of mutation cellular fraction (MCF) in small somatic mutations (SSMs) and copy number alterations (CNAs). Using a Bayesian hierarchical model, it assigns SSMs and CNAs to subclones, reconstructing tumor evolutionary trees that adhere to principles of lineage precedence, sum condition, and optional constraints based on sample presence. For deconvolution, PICTographPlus integrates tumor clonal tree structures with clone proportions across samples to resolve bulk gene expression data. It optimizes an objective function that minimizes discrepancies between observed and predicted sample-level gene expression while imposing a smoothness penalty, ensuring that closely related clones display greater gene expression similarity. Lastly, the tool conducts pathway enrichment analysis to identify statistically significant alterations in pathways connecting tumor clones.

## Installation

### JAGS
PICTographPlus uses the JAGS library for Bayesian data analysis, which is installed outside of R. JAGS can be downloaded and installed for your OS [here](https://mcmc-jags.sourceforge.io).

### PICTographPlus

```
devtools::install_github("KarchinLab/pictographPlus", build_vignettes = TRUE)
```
## Tutorial

### vignette

Detailed tutorial for using PICTographPlus is available through vignette. 

```
library(pictograph)
vignette("pictograph", package = "pictograph")
```