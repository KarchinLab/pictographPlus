#' @export
#' @param mutation_file a csv file that include information for SSMs.
#' @param rna_file bulk RNA file in integer read counts; rows are samples and columns are genes
#' @param outputDir output directory for saving all files.
#' @param copy_number_file a csv file that include information for CNA.
#' @param SNV_file a csv file that include information for germline heterozygous SNVs.
#' @param lambda weights used in deconvolution step; default is 0.2
#' @param GSEA whether to perform GSEA analysis; default is TRUE
#' @param GSEA_file geneset file in MSigDB .gmt format; the geneset name will show up in plotting
#' @param sample_presence whether to use sample presence to separate the mutations. Not applicable if dual_model is set to FALSE and a copy number file is provided.
#' @param score scoring function to estimate the number of clusters. silhouette or BIC.
#' @param max_K user defined maximum number of clusters.
#' @param min_mutation_per_cluster minumum number of mutations in each cluster.
#' @param n.iter number of iterations by JAGS.
#' @param n.burn number of burns by JAGS.
#' @param thin number of thin by JAGS.
#' @param mc.cores number of cores to use for parallel computing; not applicable to windows.
#' @param inits additional parameters by JAGS.
#' @param cluster_diff_thresh threshold to merge two clusters.
runPICTographPlus <- function(
    mutation_file,
    rna_file,
    outputDir=NULL,
    copy_number_file=NULL,
    SNV_file=NULL,
    lambda=0.2,
    GSEA = TRUE,
    GSEA_file = NULL,
    top_K = 5,
    n_permutations=1000,
    sample_presence=TRUE,
    score="silhouette", # either BIC or silhouette
    max_K = 10, 
    min_mutation_per_cluster=5, 
    cluster_diff_thresh=0.05,
    n.iter=5000, 
    n.burn=1000, 
    thin=10, 
    mc.cores=8, 
    inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
    driverFile = NULL,
    cytobandFile = NULL,
    threshes=NULL, # placeholder
    LOH = FALSE, # placeholder
    dual_model=TRUE, # placeholder
    ploidy=2, # placeholder
    pval=0.05, # placeholder
    alt_reads_thresh = 0, # placeholder
    vaf_thresh = 0, # placeholder
    cnv_max_dist=1000000, # placeholder
    cnv_max_percent=0.30, # placeholder
    tcn_normal_range=c(1.7, 2.3), # placeholder
    smooth_cnv=F, # placeholder
    autosome=T # placeholder
) {
  
  runPictograph(mutation_file,
           copy_number_file=copy_number_file,
           SNV_file=SNV_file,
           outputDir=outputDir,
           sample_presence=sample_presence,
           dual_model=dual_model, # placeholder; dual_model=FALSE still require testing
           score=score, # either BIC or silhouette
           ploidy=ploidy, # placeholder
           pval=pval, # placeholder
           max_K=max_K, 
           min_mutation_per_cluster=min_mutation_per_cluster, 
           cluster_diff_thresh=cluster_diff_thresh,
           n.iter=n.iter, 
           n.burn=n.burn, 
           thin=thin, 
           mc.cores=mc.cores, 
           inits=inits,
           threshes=threshes,
           LOH=LOH,
           driverFile=driverFile,
           cytobandFile = cytobandFile,
           alt_reads_thresh = alt_reads_thresh, # placeholder
           vaf_thresh = vaf_thresh, # placeholder
           cnv_max_dist=cnv_max_dist, # placeholder
           cnv_max_percent=cnv_max_percent, # placeholder
           tcn_normal_range=tcn_normal_range, # placeholder
           smooth_cnv=smooth_cnv, # placeholder
           autosome=autosome # placeholder
  )
  
  treeFile = paste(outputDir, "tree.csv", sep="/")
  proportionFile = paste(outputDir, "subclone_proportion.csv", sep="/")
  purityFile = paste(outputDir, "purity.csv", sep="/")
  
  X_optimal = runDeconvolution(rna_file = rna_file,
                treeFile = treeFile,
                proportionFile = proportionFile,
                purityFile = purityFile,
                outputDir = outputDir,
                lambda=lambda)
  
  
  if (GSEA) {
    # X_optimal <- read.csv(paste0(outputDir, "/clonal_expression.csv"), row.names = 1)
    runGSEA(X_optimal, outputDir, treeFile, GSEA_file=GSEA_file, top_K=top_K,
            n_permutations=n_permutations)
  }

}







