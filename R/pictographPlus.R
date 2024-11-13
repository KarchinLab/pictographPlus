#' @import GSVA, pheatmap, limma
#' @export
#' 
runPICTographPlus <- function(
    mutation_file,
    copy_number_file=NULL,
    SNV_file=NULL,
    rna_file=NULL,
    outputDir=NULL,
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
  
  mcmcMain(mutation_file,
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
  
  deconvolution(rna_file = rna_file,
                treeFile = treeFile,
                proportionFile = proportionFile,
                purityFile = purityFile,
                outputDir = outputDir)
  
  # read expression data
  X_optimal <- read.csv(paste0(outputDir, "/clonal_expression.csv"))
  rownames(X_optimal) <- X_optimal[,1]
  X_optimal <- X_optimal[,-1]
  # X_optimal <- X_optimal[1,]
  X <- mapply(as.matrix(X_optimal), FUN=as.numeric)
  X <- matrix(X, ncol=ncol(X_optimal))
  rownames(X) <- rownames(X_optimal)
  colnames(X) <- colnames(X_optimal)
  X <- t(X)
  
  # GSEA gene set
  lines <- readLines("./inst/extdata/h.all.v2024.1.Hs.symbols.gmt.txt") 
  
  gene_list <- list()
  
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    gene_set_name <- elements[1]
    genes <- elements[-c(1,2)]
    gene_list[[gene_set_name]] <- genes
  }
  
  # run GSVA
  
  gsvaPar <- gsvaParam(X, gene_list)
  gsvaPar <- ssgseaParam(X, gene_list)
  gsva.es <- gsva(gsvaPar)
  
  pheatmap(gsva.es,
           main = "Single Sample GSEA Scores",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "none",
           color = colorRampPalette(c("blue", "white", "red"))(100))
  
  
  # differetial GSEA
  data <- t(X_optimal)
  normal_sample <- data[, "0", drop = FALSE]
  disease_samples <- data[, colnames(data) != "0"]
  gene_sets <- getGmt("./inst/extdata/h.all.v2024.1.Hs.symbols.gmt.txt", geneIdType=SymbolIdentifier())

  combined_samples <- cbind(normal_sample, disease_samples)
  
  # Perform ssGSEA
  gsea_results <- gsva(gsvaParam(combined_samples, gene_sets))
  
  group <- factor(c("Normal", rep("Disease", ncol(disease_samples))))
  design <- model.matrix(~ group)
  
  # Fit the linear model
  fit <- lmFit(gsea_results, design)
  fit <- eBayes(fit)
  
  # Get the top differentially enriched gene sets
  results <- topTable(fit, coef = "groupDisease", adjust = "BH", number = Inf)
  
  
}







