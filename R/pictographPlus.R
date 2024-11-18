#' @import GSVA, pheatmap, limma
#' @export
#' 
runPICTographPlus <- function(
    mutation_file,
    rna_file,
    outputDir=NULL,
    copy_number_file=NULL,
    SNV_file=NULL,
    lambda=0.2,
    GSEA = TRUE,
    GSEA_file = NULL,
    sample_presence=TRUE,
    dual_model=TRUE, # placeholder; dual_model=FALSE still require testing
    score="silhouette", # either BIC or silhouette
    ploidy=2, # placeholder
    pval=0.05, # placeholder
    max_K = 10, 
    min_mutation_per_cluster=5, 
    cluster_diff_thresh=0.05,
    n.iter=5000, 
    n.burn=1000, 
    thin=10, 
    mc.cores=8, 
    inits=list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 123),
    threshes=NULL,
    LOH = FALSE,
    alt_reads_thresh = 0, # placeholder
    vaf_thresh = 0, # placeholder
    cnv_max_dist=2000, # placeholder
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
  
  runDeconvolution(rna_file = rna_file,
                treeFile = treeFile,
                proportionFile = proportionFile,
                purityFile = purityFile,
                outputDir = outputDir,
                lambda=lambda)
  
  
  if (GSEA) {
    X_optimal <- read.csv(paste0(outputDir, "/clonal_expression.csv"))
    runGSEA(X_optimal, outputDir, GSEA_file=GSEA_file)
  }

}







