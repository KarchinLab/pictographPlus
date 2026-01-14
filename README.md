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
- Deconvolves bulk gene expression data using **tumor clonal tree structures**.
- Performs **pathway enrichment analysis** to highlight significant transcriptomic alterations.

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
## Data and Scripts for the manuscript

All precessed data are available on Mendeley.

The scripts used for the manuscript can be found at https://github.com/KarchinLab/PICTographPlus_manuscript.
