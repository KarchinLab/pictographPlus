# pictograph2
PICTograph2 infers the clonal evolution of tumors from single or multi-region sequencing data. The tool models uncertainty of mutation cellular fraction (MCF) in small somatic mutations (SSMs) and copy number alterations (CNAs), assigning SSMs and CNAs to subclones using a Bayesian hierarchical model, and reconstruct tumor evolutionary trees that are constrained based on principles of lineage precedence, sum condition, and optionally by sample-presence. The inputs to PICTograph2 are variant ("alt allele") read counts of SSMs, sequencing depth at the SSM loci, and the somatic CNA of the tumor genome. Tumor purity and the number of mutant alleles (multiplicity) are optional parameters. If multiple tumor samples are considered, an option is available to restrict the number of possible evolutionary trees by partitioning mutations according to sample presence. PICTograph2 summarizes the posterior distributions of the mutation cluster assignments and the MCFs for each cluster by the mode. The estimates of cluster MCFs are then used to determine the most probable trees. Multiple trees that share the same score can be summarized as an ensemble tree, where edges are weighted by their concordance among constituent trees in the ensemble.

## Installation

### JAGS
PICTograph2 uses the JAGS library for Bayesian data analysis, which is installed outside of R. JAGS can be downloaded and installed for your OS [here](https://mcmc-jags.sourceforge.io).

### PICTograph2

```
devtools::install_github("KarchinLab/pictograph2", build_vignettes = TRUE)
```
## Tutorial

### vignette

A vignette that includes a toy example is available in R. It is similar to the Tutorial listed below.

```
library(pictograph2)
vignette("pictograph2", package = "pictograph2")
```

### parameters for each function
The details of the parameters for each function can be viewed in R using help(function). For example:
```
help(mcmcMain)
```

### 1. Input data

Pictograph2 takes input data in multiple formats for flexible user inputs:

1) A single csv file that contains SSM and CNA information

* The first option is to provide a single csv file that contains at least columns named "sample", "mutation", "total_reads", "alt_reads", "tumor_integer_copy_number", and "cncf". Example input files can be found under "inst/extdata/examples". See files that starts with example1.

    | sample | mutation | total_reads | alt_reads | tumor_integer_copy_number | cncf |
    | ---- | ---- | ---- | ---- | ---- | ---- |
    | sample1 | mut1 | 100 | 67 | 4 | 0.8 |
    | sample1 | mut2 | 50 | 67 | 3 | 1 |
    | sample2 | mut1 | 100 | 50 | 4 | 0.7 |

* Users can also provide an optional column "major_integer_copy_number" that provides the information of the integer copy number of the major allele. If "major_integer_copy_number" is not provided, it will be estimated using an internal function built in the package. Example input files can be found under "inst/extdata/examples". See files that starts with example2.

    | sample | mutation | total_reads | alt_reads | tumor_integer_copy_number | major_integer_copy_number | cncf |
    | ---- | ---- | ---- | ---- | ---- | ---- | --- |
    | sample1 | mut1 | 100 | 67 | 4 | 3 | 0.8 |
    | sample1 | mut2 | 50 | 67 | 3 | 2| 1 |
    | sample2 | mut1 | 100 | 50 | 4 | 3| 0.7 |

* Another optional column is "purity" column that provides the information of normal contamination of a sample. Putiry of 0.8 wil be used if not provided. 

    | sample | mutation | total_reads | alt_reads | tumor_integer_copy_number | major_integer_copy_number | cncf | purity |
    | ---- | ---- | ---- | ---- | ---- | ---- | --- | --- |
    | sample1 | mut1 | 100 | 67 | 4 | 3 | 0.8 | 0.9 |
    | sample1 | mut2 | 50 | 67 | 3 | 2| 1 | 0.9 |
    | sample2 | mut1 | 100 | 50 | 4 | 3| 0.7 | 0.9 |

2) Two csv files, one for SSM and one for CNA

* The second option is to provide the SSM read counts and copy number alterations (CNA) in two separate files. In this case, the SSM file should contain columns "sample", "mutation", "total_reads", "alt_reads", "chrom", "start", and "end". The "purity" column with be optional. The CNA file should contain columns "sample", "chrom", "start", "end", "tcn" and "baf". See files that starts with example3.

3) Three csv files, one for ssm, one for CNA, and one for germline heterozygous single nucleotide variants (SNV)

* In the absence of the baf column, users should provide an additional file that contains the count of the heterozygous germline SNVs that has the information about the "chroms", the "position" of the germline heterozygous SNV, "ref" and "alt" allele, and the reference and altenative reads counts in germline (normal) sample as well as all other samples. Note: the sample name should matched the sample name used in the SSM and CNA file. See files that starts with example4.

### 2. Run PICTograph2 using the main function

Run an automated pipeline of the tool using the function mcmcMain. The only required user input is the "mutation_file", along with the "copy_number_file" and "SNV_file" if using format 2 or 3 mentioned above. Users can specify the destination of output files using "outputDir" option. If the outputDir is not specified, output files will be stored in the current working directory.

```
# here we use the example1 file as a demo
mcmcMain(mutation_file=system.file('extdata/examples/example1_snv.csv', package = 'pictograph2'), outputDir=system.file('extdata/examples/output', package = 'pictograph2'))
```

This will run PICTograph2 and save the output files in 'extdata/examples/output' directory. 

### 3. Parameters

| Parameter | Description | type | options |
| ------------- | ------------- | ------- | ------ |
| mutation_file | a csv file that include information for SSMs. Details explained in the Input data section. See inst/extdata/examples for an example. | string | Required 
| copy_number_file | a csv file that include information for CNA. Details explained in the Input data section. See inst/extdata/examples for an example. | string | default: NULL
| SNV_file | a csv file that include information for germline heterozygous SNVs. Details explained in the Input data section. See inst/extdata/examples for an example. | string | default: NULL
| outputDir | output directory for saving all files. default: current directory. | string | default: NULL
| dual_model | whether to use one model or two separate models. Details can be found in Models section. Applicable if a copy number file is provided. | boolean| default: TRUE
| sample_presence | whether to use sample presence to separate the mutations. Not applicable if dual_model is set to FALSE and a copy number file is provided. | boolean | default: TRUE
| score | scoring function to estimate the number of clusters. | string | silhouette or BIC. default: silhouette
| max_K | user defined maximum number of clusters. | integer | default: 10
| min_mutation_per_cluster | minumum number of mutations in each cluster. | integer | default: 5
| cluster_diff_thresh | threshold to merge two clusters. | float | default: 0.05
| n.iter | number of iterations by JAGS.| integer | default: 5000
| n.burn | number of burns by JAGS. | integer | default: 1000
| thin | number of thin by JAGS. | integer | default: 10
| inits | additional parameters by JAGS. | list | default: list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123)
| mc.cores | number of cores to use for parallel computing; not applicable to windows. | integer | default: 8


### 4. Output files

| File name | Description |
| ------------- | ------------- |
| mcf.csv | estimated MCF for each cluster in each sample  |
| clusterAssign.csv | cluster assignment of each SSM/CNA to each cluster  |
| CN_results.csv | estimation of the integer and major copy number of CNAs; will be empty if the input file is a single csv file (i.e. no separate CNA information)|
| tree.csv | the tree with the highest score; all trees with tied highest score is available under all_trees directory. |
| subclone_proportion.csv | estimated proportion of each cluster in each sample |
| purity.csv | estimated purity for each sample, based on the tree structure |
| tree_ensemble.png | the image of the ensemble tree of all trees with the same highest score  |
| tree.png | the image of a tree with the best score  |
| upsetR.png | the mutation profiles between samples; only available if number of samples is bigger than 1.|
| violin.png | a violin plot of the MCF|
| mcf.png | the MCF chain trace from JAGS|

### 5. Models