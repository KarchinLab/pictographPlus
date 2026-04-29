#' run both tumor evolution reconstruction and clone-specific transcriptomic profile deconvolution at once
#' 
#' Tumor evolution drives intra-tumor heterogeneity, yet functional clone-specific profiles remain poorly understood. 
#' We developed PICTographPlus, a computational tool that integrates bulk DNA and RNA sequencing to reconstruct clone-specific transcriptomic profiles 
#' and uncover transcriptional transitions between clones, offering a powerful framework for cancer evolution studies.
#'  
#' @export
#' 
#' @param mutation_file a csv file that include information for SSMs. See vignette for details.
#' @param rna_file bulk RNA file in integer read counts; rows are genes and columns are samples (first column = gene IDs). See vignette for details.
#' @param copy_number_file a csv file that include information for CNA. See vignette for details.
#' @param SNV_file a csv file that include information for germline heterozygous SNVs. See vignette for details.
#' @param outputDir output directory for saving all files.
#' @param sample_presence whether or not to use sample presence to separate the mutations; default: TRUE 
#' @param score scoring function to estimate the number of clusters. silhouette or BIC; default: silhuette
#' @param max_K user defined maximum number of clusters; default: 10
#' @param min_mutation_per_cluster minumum number of mutations in each cluster; default: 5
#' @param min_cluster_thresh minimum MCF for each cluster; default: 0.05
#' @param cluster_diff_thresh difference threshold to merge two clusters: default: 0.05
#' @param n.iter number of iterations by JAGS; default: 5000
#' @param n.burn number of burns by JAGS; default: 1000
#' @param thin number of thin by JAGS; default: 10
#' @param mc.cores number of cores to use for parallel computing; not applicable to windows; default: 8
#' @param inits additional parameters by JAGS.
#' @param LOH whether or not to include copy number segments that are copy neutral but LOH; default: FALSE
#' @param purity_min minimum purity for tumor samples; default: 0.2
#' @param driverFile list of driver genes used for visualization. See vignette for details.
#' @param cytobandFile list of cytoband regions used for visualization. See vignette for details.
#' @param alt_reads_thresh minimum number of alternative read count for a SSM to be included in the analysis; default: 0
#' @param vaf_thresh minimum VAF for a SSM to be included in the analysis; default: 0
#' @param tcn_normal_range range of total copy number considered as copy-neutral; default: c(1.75,2.3)
#' @param filter_cnv whether or not to filter copy number alterations; default: TRUE
#' @param smooth_cnv whether or not to process copy number alterations across samples to unify the segment start and end postions; default: TRUE
#' @param autosome to only include autosomes; default: TRUE
#' @param cnv_min_length minimum length of copy number alterations for it to be included in analysis
#' @param lambda regularisation strength for deconvolution; default: 0.01 (elastic_net star-best)
#' @param use_star_tree if TRUE (default), deconvolve with a star topology (root → all clones)
#'   rather than the inferred tree. Recommended until topology advantage is confirmed.
#' @param model deconvolution model; one of "elastic_net" (default), "tree_delta", "plain",
#'   "adaptive", "adaptive_v2", "plain_debiased", "fused_ew". See \code{runDeconvolution}.
#' @param GSEA whether to perform GSEA analysis; default is TRUE
#' @param GSEA_file geneset file in MSigDB .gmt format; the geneset name will show up in plotting
#' @param top_K top_K significant pathways to be plotted as GSEA results; default: 5
#' @param n_permutations number of permutations in fgsea; default: 10000
#' @param purityFile purity file for each sample; default: NULL, in which case the
#'   \code{purity.csv} written by \code{runPictograph} in \code{outputDir} is used
#'   automatically if present. Supply a path to override.
#' @param normalize normalize the raw count using DESeq2; default: TRUE
runPICTographPlus <- function(
    mutation_file,
    rna_file,
    outputDir=NULL,
    copy_number_file=NULL,
    SNV_file=NULL,
    lambda=0.01,
    use_star_tree=TRUE,
    model="elastic_net",
    GSEA = TRUE,
    GSEA_file = NULL,
    top_K = 5,
    normalize=TRUE,
    purityFile = NULL,
    n_permutations=10000,
    sample_presence=FALSE,
    score="BIC", # either BIC or silhouette
    max_K = 10, 
    min_mutation_per_cluster=1, 
    min_cluster_thresh=0.05, 
    cluster_diff_thresh=0.05,
    n.iter=5000, 
    n.burn=1000, 
    thin=10, 
    mc.cores=8, 
    inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
    driverFile = NULL,
    cytobandFile = NULL,
    LOH = FALSE,
    purity_min=0.2, 
    alt_reads_thresh = 0, 
    vaf_thresh = 0, 
    cnv_min_length=1000000,
    tcn_normal_range=c(1.75, 2.3),
    filter_cnv = T, 
    smooth_cnv = T, 
    autosome = T,
    threshes = NULL,
    dual_model = TRUE, 
    pval = 0.05, 
    ploidy = 2
) {
  
  runPictograph(mutation_file,
           copy_number_file=copy_number_file,
           SNV_file=SNV_file,
           outputDir=outputDir,
           sample_presence=sample_presence,
           dual_model=dual_model,
           score=score,
           ploidy=ploidy,
           pval=pval,
           max_K=max_K, 
           min_mutation_per_cluster=min_mutation_per_cluster, 
           min_cluster_thresh=min_cluster_thresh,
           cluster_diff_thresh=cluster_diff_thresh,
           n.iter=n.iter, 
           n.burn=n.burn, 
           thin=thin, 
           mc.cores=mc.cores, 
           inits=inits,
           threshes=threshes,
           LOH=LOH,
           purity_min=purity_min,
           selectedMutFile = driverFile,
           cytobandFile = cytobandFile,
           alt_reads_thresh = alt_reads_thresh,
           vaf_thresh = vaf_thresh,
           cnv_min_length = cnv_min_length, 
           tcn_normal_range=tcn_normal_range,
           filter_cnv=filter_cnv,
           smooth_cnv=smooth_cnv, 
           autosome=autosome
  )
  
  treeFile = paste(outputDir, "tree.csv", sep="/")
  proportionFile = paste(outputDir, "subclone_proportion.csv", sep="/")

  # Default: use the purity.csv written by runPictograph (MCMC-main.R).
  # Users can still override by passing purityFile explicitly.
  if (is.null(purityFile)) {
    auto_purity <- file.path(outputDir, "purity.csv")
    if (file.exists(auto_purity)) {
      purityFile <- auto_purity
      message("Using purity file written by runPictograph: ", auto_purity)
    }
  }

  X_optimal = runDeconvolution(rna_file       = rna_file,
                treeFile       = if (use_star_tree) NULL else treeFile,
                proportionFile = proportionFile,
                normalize      = normalize,
                purityFile     = purityFile,
                outputDir      = outputDir,
                lambda         = lambda,
                use_star_tree  = use_star_tree,
                model          = model)

  if (GSEA) {
    # GSEA uses the DNA-derived tree regardless of whether deconvolution used a star tree.
    # Edge-wise enrichment must follow the biological clonal phylogeny.
    runGSEA(X_optimal, outputDir, treeFile, GSEA_file=GSEA_file, top_K=top_K,
            n_permutations=n_permutations)
  }

}







