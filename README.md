
# **PICTographPlus**

PICTographPlus is a computational tool that integrates bulk DNA and RNA sequencing to reconstruct clone-specific transcriptomic profiles and uncover transcriptional transitions between tumor clones, providing a robust framework for studying cancer evolution. It comprises three core models: 1) tumor evolution inference using genomic data (runPictograph), 2) bulk RNA expression deconvolution based on tumor evolution (runDeconvolution), and 3) gene set enrichment analysis (GSEA) to examine transcriptomic differences between clones (runGSEA).

The tool infers tumor clonal evolution from single or multi-region sequencing data by modeling the uncertainty of mutation cellular fraction (MCF) in small somatic mutations (SSMs) and copy number alterations (CNAs). Using a Bayesian hierarchical model, it assigns SSMs and CNAs to subclones, reconstructing tumor evolutionary trees that adhere to principles of lineage precedence, sum condition, and optional constraints based on sample presence. For deconvolution, PICTographPlus integrates tumor clonal tree structures with clone proportions across samples to resolve bulk gene expression data. It optimizes an objective function that minimizes discrepancies between observed and predicted sample-level gene expression while imposing a smoothness penalty, ensuring that closely related clones display greater gene expression similarity. Lastly, the tool conducts pathway enrichment analysis to identify statistically significant alterations in pathways connecting tumor clones.


## **Overview**

PICTographPlus is a computational tool that integrates bulk DNA and RNA sequencing data to:

1. **Reconstruct Clone-Specific Transcriptomic Profiles**
2. **Infer Tumor Evolution**
3. **Identify Transcriptional Transitions Between Clones**

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

--- 

## **Access Tutorial**
```r
library(pictographPlus)
vignette("pictographPlus", package = "pictographPlus")
```

---

## Input data for genomic data for tumor evolution reconstruction

PICTographPlus takes input data in multiple formats for flexible user inputs:

1. Three CSV files, one each for SSMs, CNAs, and germline heterozygous SNVs (RECOMMENDED)
2. Two csv files, one for SSM and one for CNA
3. A single csv file that contains SSM and CNA information. 

### 1) Three CSV files, one each for SSMs, CNAs, and germline heterozygous SNVs (RECOMMENDED)

The recommended input for tumor evolution reconstruction is to provide individual files for SSMs, CNAs and germline heterozygous SNVs. The SSM file should contain at least sample, mutation, total_reads, alt_reads, chrom, start, and end columns, with the purity column being optional. The CNA file should contain sample, chrom, start, end, and tcn columns. The total copy number (tcn) can be inferred from many copy number callers. In case a copy number caller outputs the $log2$ ratio of a segment, the tcn can be calculate using $tcn = 2^{\log_2 R}$. 

The SNV file contains the count of the heterozygous germline SNVs that has the information about the "chrom", the "position" of the germline heterozygous SNV, "ref" and "alt" allele, and the reference and altenative reads counts in germline (normal) sample as well as all other samples. Note: the sample name (sample1, sample2, ... etc.) should matched the sample name used in the SSM and CNA file. The germline positions can be obtained using tools such as ```GATK HaplotypeCaller```. The reference an alternative read counts for each tumor sample can be obtained using ```samtools mpileup```. We provide a python script, ```getPileUp.py```, to get the pileup information. 

```python
python getPileUp.py -v germline.vcf -b sample1.bam sample2.bam ... -o outputDir -f hg38.fa [--submit]
```
| parameter | description | option |
| - | - | - |
| -v| vcf output from HaplotypeCaller | required
| -b| tumor bam files for tumor samples | required
| -o| output directory | required
| -f| human reference genome | required
| --submit | submit as a slurm job | optional

Examples of the three csv files are:

* SSM file

    | sample | mutation | total_reads | alt_reads | chrom | start | end | purity (optional) |
    | ---- | ---- | ---- | ---- | ---- | ---- | --- | --- |
    | sample1 | mut1 | 100 | 67 | chr1 | 1000 | 1000 | 0.9
    | sample1 | mut2 | 50 | 67 | chr2 | 3000| 3000 | 0.9
    | sample2 | mut1 | 100 | 50 | chr1 | 1000 | 1000 | 0.9

* CNA file
    | sample | chrom | start | end | tcn | 
    | ---- | ---- | ---- | ---- | ---- |
    | sample1 | chr1 | 10 | 10000 | 3.6 |
    | sample1 | chr2 | 2000| 3000 | 3.4 |
    | sample2 | chr1 | 10 | 10000 | 3.6 |

* SNV file of germline heterozygous SNVs

    | chrom | position | ref | alt | germline_ref | germline_alt | sample1_ref | sample1_alt | sample2_ref | sample2_alt |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | chr1 | 20 | A | C | 50 | 50 | 70 | 30 | 60 | 40 |
    | chr1 | 50 | G | T | 55 | 45 | 76 | 34 | 80 | 20 |

### 2) Two csv files, one for SSM and one for CNA

The second option is to provide the SSM read counts and copy number alterations (CNA) in two separate files. In this case, the SSM file should contain columns "sample", "mutation", "total_reads", "alt_reads", "chrom", "start", and "end". The "purity" column with be optional. The CNA file should contain columns "sample", "chrom", "start", "end", "tcn" and "baf". 

* SSM file

    | sample | mutation | total_reads | alt_reads | chrom | start | end |
    | ---- | ---- | ---- | ---- | ---- | ---- | --- |
    | sample1 | mut1 | 100 | 67 | chr1 | 1000 | 1000 |
    | sample1 | mut2 | 50 | 67 | chr2 | 3000| 3000 |
    | sample2 | mut1 | 100 | 50 | chr1 | 1000 | 1000 |

* CNA file with baf column

    | sample | chrom | start | end | tcn | baf |
    | ---- | ---- | ---- | ---- | ---- | ---- |
    | sample1 | chr1 | 10 | 10000 | 3.6 | 0.3 |
    | sample1 | chr2 | 2000| 3000 | 3.4 | 0.3 |
    | sample2 | chr1 | 10 | 10000 | 3.6 | 0.4 |


### 3) A single csv file that contains SSM and CNA information. 

The last option is to provide a single csv file that contains at least columns named "sample", "mutation", "total_reads", "alt_reads", "tumor_integer_copy_number", and "cncf". Set cncf to 0 if a mutation has no copy number alteration. Users can also provide an optional column "major_integer_copy_number" that provides the information of the integer copy number of the major allele. If "major_integer_copy_number" is not provided, it will be estimated using an internal function built in the package. Another optional column is "purity" column that provides the information of normal contamination of a sample. 

NOTE: using this option will generate trees with SSMs only, CNA will not be assigned to clusters but only used for VAF correction.

* SSM file with minimal information

    | sample | mutation | total_reads | alt_reads | tumor_integer_copy_number | cncf | major_integer_copy_number (optional) | purity (optional)
    | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
    | sample1 | mut1 | 100 | 67 | 4 | 0.8 | 3 | 0.8
    | sample1 | mut2 | 100 | 67 | 4 | 0.8 | 2 | 0.8
    | sample2 | mut1 | 100 | 40 | 2 | 1 | 3 | 0.8

---

## Input data for bulk RNA expression

The RNA file should be a csv file of columns Gene, followed by the tumor samples (tumor sample name should match that of the genomic input), and lastly the read counts of a matched normal sample. The read counts should be normalized to transcript per million (TPM).

| Gene | sample1 | sample2 | sampleN 
| ---- | ---- | ---- | ---- |
| gene1 | 1 | 1| 3
| gene2 | 10 | 12 | 5
| gene3 | 1 | 0 | 0

---

## Running PICTographPlus in one step

PICTographPlus can be run using the function ```runPICTographPlus```, which runs both tumor evolution reconstruction and clone-specific transcriptomic profile deconvolution. The required files include files for genomic data and RNA expression data. 

```
runPICTographPlus(mutation_file, copy_number_file, SNV_file, rna_file, outputDir)
```

Detailed documentation for the function can be found:
```
help(runPICTographPlus)
```

### Parameters for ```runPICTographPlus```

| Parameter | Description | default |
| ------------- | ------------- | ------ |
| mutation_file | a csv file that include information for SSMs | Required
| rna_file | bulk RNA file in integer read counts; rows are samples and columns are genes | Required
| copy_number_file | a csv file that include information for CNA | NULL
| SNV_file | a csv file that include information for germline heterozygous SNVs | NULL
| outputDir | output directory for saving all files | NULL
| sample_presence | whether or not to use sample presence to separate the mutations | TRUE 
| score | scoring function to estimate the number of clusters. silhouette or BIC | silhuette
| max_K | user defined maximum number of clusters | 10
| min_mutation_per_cluster | minumum number of mutations in each cluster | 5
| min_cluster_thresh | minimum MCF for each cluster | 0.05
| cluster_diff_thresh | difference threshold to merge two clusters | 0.05
| n.iter | number of iterations by JAGS | 5000
| n.burn | number of burns by JAGS | 1000
| thin | number of thin by JAGS | 10
| mc.cores | number of cores to use for parallel computing; not applicable to windows | 8
| inits | additional parameters by JAGS | list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123)
| LOH | whether or not to include copy number segments that are copy neutral but LOH | FALSE
| purity_min | minimum purity for tumor samples | 0.2
| driverFile | list of driver genes used for visualization | NULL
| cytobandFile | list of cytoband regions used for visualization | NULL
| alt_reads_thresh | minimum number of alternative read count for a SSM to be included in the analysis | 0
| vaf_thresh | minimum VAF for a SSM to be included in the analysis | 0
| tcn_normal_range | range of total copy number considered as copy-neutral | c(1.75,2.3)
| filter_cnv | whether or not to filter copy number alterations | TRUE
| smooth_cnv | whether or not to process copy number alterations across samples to unify the segment start and end postions | TRUE
| autosome | to only include autosomes | TRUE
| cnv_min_length | minimum length of copy number alterations for it to be included in analysis | 1000000
| lambda weights | used in deconvolution step | 0.2
| GSEA | whether to perform GSEA analysis | TRUE
| GSEA_file | geneset file in MSigDB .gmt format; the geneset name will show up in plotting | NULL
| top_K | top_K significant pathways to be plotted as GSEA results | 5
| n_permutations | number of permutations in fgsea | 10000

Output files will be stored in the outputDir a user specifies, or in the current directory if the outputDir is not provided. 

### Output files

| File name | Description |
| ------------- | ------------- |
| mcf.csv | estimated MCF for each cluster in each sample  |
| clusterAssign.csv | cluster assignment of each SSM/CNA to each cluster  |
| CN_results.csv | estimation of the integer and major copy number of CNAs |
| tree.csv | the tree with the highest score; all trees with tied highest score is available under all_trees directory. |
| subclone_proportion.csv | estimated proportion of each cluster in each sample |
| purity.csv | estimated purity for each sample, based on the tree structure |
| tree.png | the image of a tree with the best score  |
| upsetR.png | the mutation profiles between samples; only available if number of samples is bigger than 1.|
| mcf.png | the MCF chain trace from JAGS|
| mutationClusterAssign.csv | table of all mutations information in all samples
| clone_expression.csv | clone level gene expression for each clone
| GSEA | directory contains all files from GSEA analysis

---

## Running PICTographPlus in multiple steps

Alternatively, the user can run tumor evolurion reconstruction and bulk RNA deconvolution in separate steps. 

**Tumor evolution reconstruction** can be run using:
```
runPictograph(mutation_file, copy_number_file, SNV_file, outputDir)
```

**Bulk RNA exression deconvolution** can be run using:
```
runDeconvolution(rna_file, treeFile, proportionFile, purityFile, outputDir)
```
where the treeFile, proportionFile, and purityFile are outputs of ```runPictograph``` function. Users may also choose other tools to get these information.

**GSEA analysis using fgsea** can be run using:
```
X_optimal <- read.csv(paste0(outputDir, "/clonal_expression.csv), row.names=1, check.names=FALSE)
runGSEA(X_optimal, outputDir, treeFile, GSEA_file)
```
where GSEA_file is a text file of pathways of interest that can be obtained from resources such as [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp).

---

## Additional file format

**treeFile** - csv file of tumor evolution tree used by ```runDeconvolution``` and ```runGSEA```

| edge | parent | child |
| - | - | - |
| root->1 | root | 1|
| 1->2 | 1 | 2

**proportionFile** - csv file of subclone proportions used by ```runDeconvolution```. The row names are clones.

| | sample1 | sample2 |
| - | - | - |
| 1 | 0.12 | 0.3 |
| 2 | 0.82 | 0.7 |

**purityFile** - csv file of tumor purity used by ```runDeconvolution```. 

| sample1 | sample2 |
| - | - |
| 0.7 | 0.5 |

**driverFile** - csv file of with driver mutation information that will be used for plotting by ```runPictograph```. If using this file, make sure the mutation in mutation_file follows the format of $gene_$extra_information where the $gene should be the gene name. Gene_type should be either oncogene or tumor_supressor, leave blank if unknown. 

| gene | chrom | start | end | gene_type |
| - | - | - | - | - |
| KRAS | chr12 | 25205246 | 25250936 | oncogene |
| SMAD4 | chr18 | 51028528 | 51085045 | tumor_suppressor |

**cytobandFile** - tsv file of cytoband information used for plotting by ```runPictograph```.

| | | | | |
| - | - | - | - | - |
| chr1 | 0 | 2300000 | p36.33 | gneg |
| chr1 2300000 | 5300000 | p36.32 | gpos25 |


