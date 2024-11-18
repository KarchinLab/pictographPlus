#' @import igraph
#' @export
#' 
runDeconvolution <- function(rna_file,
                          treeFile,
                          proportionFile,
                          purityFile,
                          outputDir,
                          lambda=0.2){
  
  # rnaFile <- "~/Karchin Lab Dropbox/Jillian Lai/U01-Wood-Fertig/analysis/htan-mcl-pre-cancer-pancreas/semaan_analysis/pictographPlus/MCL111_001_expression.csv"
  # treeFile <- "~/Karchin Lab Dropbox/Jillian Lai/U01-Wood-Fertig/analysis/htan-mcl-pre-cancer-pancreas/pictograph2_latest/MCL111_001/tree.csv"
  # proportionFile <- "~/Karchin Lab Dropbox/Jillian Lai/U01-Wood-Fertig/analysis/htan-mcl-pre-cancer-pancreas/pictograph2_latest/MCL111_001/subclone_proportion.csv"
  # purityFile <- "~/Karchin Lab Dropbox/Jillian Lai/U01-Wood-Fertig/analysis/htan-mcl-pre-cancer-pancreas/pictograph2_latest/MCL111_001/purity.csv"
  # outputDir <- "~/Karchin Lab Dropbox/Jillian Lai/U01-Wood-Fertig/analysis/htan-mcl-pre-cancer-pancreas/semaan_analysis/pictographPlus"
  
  rnaData <- read.csv(rna_file, row.names=1)
  propData <-read.csv(proportionFile, row.names=1)
  purityData <- read.csv(purityFile)
  tree <- read.csv(treeFile)
  
  proportionDF <- propData * as.numeric(purityData[1,])
  proportionDF <- rbind(proportionDF, 1 - purityData)
  rownames(proportionDF)[nrow(proportionDF)] <- "0"
  
  lesion <- setdiff(colnames(rnaData), colnames(proportionDF))
  proportionDF[[lesion]] <- 0
  proportionDF[nrow(proportionDF), lesion] <- 1
  proportionDF <- round(proportionDF, 4)
  proportionDF <- proportionDF[, colnames(rnaData)]
  proportionDF <- proportionDF[order(rownames(proportionDF)), ]
  
  edges <- list()
  con <- file(treeFile, "r")
  on.exit(close(con))
  readLines(con, n = 1)
  
  while(TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break  # Exit loop if no more lines
    line <- strsplit(line, ",")[[1]]
    
    # Add edges based on conditions
    if (line[2] == "root") {
      edges <- append(edges, list(c(0, as.integer(line[3]))))
    } else if (line[3] == "root") {
      edges <- append(edges, list(c(as.integer(line[2]), 0)))
    } else {
      edges <- append(edges, list(c(as.integer(line[2]), as.integer(line[3]))))
    }
  }
  
  # Create a graph and add edges
  edge_list <- as.matrix(do.call(rbind, edges))  # Convert to matrix if not already
  edge_list <- apply(edge_list, 2, as.integer) 
  edge_list <- edge_list + 1
  if (is.null(ncol(edge_list))) {
    L <- matrix(0, nrow = 2, ncol = 2)
  } else {
    G <- graph_from_edgelist(edge_list, directed = FALSE)
    laplacian_mat <- laplacian_matrix(G)
    laplacian_mat[1,] <- 0
    laplacian_mat[,1] <- 0
    L <- as.matrix(laplacian_mat)
  }
  
  Y <- as.matrix(t(rnaData))  # Transpose the data frame
  pi <- as.matrix(t(proportionDF))
  
  X_optimal <- optimize_X(Y, pi, L, lambda)
  write.csv(X_optimal, file = paste0(outputDir, "/clonal_expression.csv"))
  
  # if (GSEA) {
  #   GSEA_analysis(X_optimal, outputDir, GSEA_file=GSEA_file)
  # }
  
}


#' @export
#' @import GSVA, pheatmap, limma, GSEABase
runGSEA <- function(X_optimal,
                    outputDir,
                    GSEA_file=NULL) {
  # read expression data
  # X_optimal <- read.csv(paste0(outputDir, "/clonal_expression.csv"))
  rownames(X_optimal) <- X_optimal[,1]
  X_optimal <- X_optimal[,-1]
  # X_optimal <- X_optimal[1,]
  X <- mapply(as.matrix(X_optimal), FUN=as.numeric)
  X <- matrix(X, ncol=ncol(X_optimal))
  rownames(X) <- rownames(X_optimal)
  colnames(X) <- colnames(X_optimal)
  X <- t(X)
  
  # GSEA gene set
  if (is.null(GSEA_file)) {
    GSEA_file = system.file("extdata", "h.all.v2024.1.Hs.symbols.gmt.txt", package="pictographPlus")
  }
  
  lines <- readLines(GSEA_file) 
  
  gene_list <- list()
  
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    gene_set_name <- elements[1]
    genes <- elements[-c(1,2)]
    gene_list[[gene_set_name]] <- genes
  }
  
  # run GSVA
  # gsvaPar <- gsvaParam(X, gene_list)
  gsvaPar <- ssgseaParam(X, gene_list)
  gsva.es <- gsva(gsvaPar)
  
  write.csv(gsva.es, file = paste0(outputDir, "/gsea_results.csv"))
  
  png(paste0(outputDir, "/gsea_heatmap.png"), width = 1200, height = 1000, res = 150)
  pheatmap(gsva.es,
           main = "Single Sample GSEA Scores",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "none",
           color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
  
  # # differetial GSEA
  # data <- t(X_optimal)
  # normal_sample <- data[, "0", drop = FALSE]
  # disease_samples <- data[, colnames(data) != "0"]
  # gene_sets <- getGmt(GSEA_file, geneIdType=SymbolIdentifier())
  # 
  # combined_samples <- cbind(normal_sample, disease_samples)
  # 
  # # Perform ssGSEA
  # gsea_results <- gsva(gsvaParam(combined_samples, gene_sets))
  # 
  # group <- factor(c("Normal", rep("Disease", ncol(disease_samples))))
  # design <- model.matrix(~ group)
  # 
  # # Fit the linear model
  # fit <- lmFit(gsea_results, design)
  # fit <- eBayes(fit)
  # 
  # # Get the top differentially enriched gene sets
  # results <- topTable(fit, coef = "groupDisease", adjust = "BH", number = Inf)
  
}

optimize_X <- function(Y, pi, L, lambda_=0.1, X_init = NULL, learning_rate = 0.01, max_iter = 100000, tol = 1e-10, round_result = TRUE) {
  # Get dimensions of pi and Y
  k <- ncol(pi)
  n <- ncol(Y)
  
  # Initialize X
  if (is.null(X_init)) {
    X <- matrix(rnorm(k * n), nrow = k, ncol = n)
  } else {
    X <- pmax(X_init, 0)
  }
  
  X_new <- X  # Copy X for the iteration
  
  for (i in 1:max_iter) {
    # Gradient calculation
    gradient <- -2 * t(pi) %*% (Y - pi %*% X_new) + 2 * lambda_ * L %*% X_new
    
    # Gradient descent step
    X_new <- X_new - learning_rate * gradient
    
    # Projection step: Enforce non-negativity
    X_new <- pmax(X_new, 0)
    
    # Check for convergence
    diff <- norm(X_new - X, "F")  # Frobenius norm
    if (diff < tol) {
      cat("Converged after", i, "iterations.\n")
      break
    }
    
    # Optional logging every 10000 iterations
    if (i %% 10000 == 0) {
      cat("diff is", log10(diff), "after iteration", i, "\n")
    }
    
    X <- X_new  # Update X to the new value for the next iteration
  }
  
  # Return rounded or raw result based on round_result flag
  if (!round_result) {
    return(X_new)
  }
  return(round(X_new))
  
}
