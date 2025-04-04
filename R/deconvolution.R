#' bulk RNA deconvolution using tumor evolution information
#' 
#' @export
runDeconvolution <- function(rna_file,
                             treeFile,
                             proportionFile,
                             outputDir,
                             purityFile=NULL,
                             lambda=0.2){
  rnaData <- read.csv(rna_file, row.names=1)
  propData <-read.csv(proportionFile, row.names=1)
  tree <- read.csv(treeFile)
  
  if (is.null(purityFile)) {
    proportionDF <- propData
    proportionDF <- rbind(proportionDF, 0)
  } else {
    purityData <- read.csv(purityFile)
    proportionDF <- t(t(propData) * as.numeric(purityData[1, ]))
    proportionDF <- rbind(proportionDF, 1 - purityData)
  }
  
  rownames(proportionDF)[nrow(proportionDF)] <- "0"
  
  lesion <- setdiff(colnames(rnaData), colnames(proportionDF))
  
  if (length(lesion) > 0) {
    proportionDF[[lesion]] <- 0
    proportionDF[nrow(proportionDF), lesion] <- 1
  }
  
  proportionDF <- round(proportionDF, 4)
  proportionDF <- proportionDF[, colnames(rnaData)]
  proportionDF <- proportionDF[order(rownames(proportionDF)), ]
  
  edge_list <- read_tree(treeFile)
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
  
  return(X_optimal)
}



#' GSEA analysis using fgsea
#' 
#' @export
#' @import ggplot2 fgsea ggrepel pheatmap DESeq2
runGSEA <- function(X_optimal,
                    outputDir,
                    treeFile,
                    GSEA_file=NULL,
                    top_K=5,
                    n_permutations=1000) {
  
  GSEA_dir = paste0(outputDir, "/GSEA")
  suppressWarnings(dir.create(GSEA_dir))
  
  X <- mapply(as.matrix(X_optimal), FUN=as.numeric)
  X <- matrix(X, ncol=ncol(X_optimal))
  rownames(X) <- rownames(X_optimal)
  colnames(X) <- colnames(X_optimal)
  X <- t(X)
  
  # GSEA gene set
  if (is.null(GSEA_file)) {
    GSEA_file = system.file("extdata", "h.all.v2024.1.Hs.symbols.gmt.txt", package="pictographPlus")
  }
  
  gene_list <- read_GSEA_file(GSEA_file)
  
  edge_list <- read_tree(treeFile)
  
  gsea_results_list <- list()
  
  if (is.null(ncol(edge_list))) {
    
    edge_list <- as.character(edge_list)
    
    sample1 = edge_list[1]
    sample2 = edge_list[2]
    
    gsea_results_list[[1]] <- GSEA_diff(X, sample1, sample2, gene_list, GSEA_dir, n_permutations, top_K)
    
    leadingEdgePlot(log2(X+1), sample1, sample2, GSEA_dir)
    
  } else {
    edge_list <- apply(edge_list, 2, as.character)
    
    for (i in 1:nrow(edge_list)) {
      
      sample1 = edge_list[i, 1]
      sample2 = edge_list[i, 2]
      
      gsea_results_list[[i]] <- GSEA_diff(X, sample1, sample2, gene_list, GSEA_dir, n_permutations, top_K)

      leadingEdgePlot(log2(X+1), sample1, sample2, GSEA_dir)
    }
  }
  
}

#'@import purrr
leadingEdgePlot <- function(X, sample1, sample2, GSEA_dir) {
  
  pdf_filename <- paste0(GSEA_dir, "/", "heatmaps_", sample1, "_", sample2, "_top30.pdf")
  pdf(pdf_filename, width = 8, height = 6)  # Open PDF device
  
  pathway_data <- read.csv(paste0(GSEA_dir, "/", "clone", sample2, "_", sample1, "_GSEA_diff.csv"))
  
  filtered_pathways <- pathway_data %>%
    filter(padj < 0.05) %>%
    mutate(leadingEdge = strsplit(as.character(leadingEdge), ";")) %>%
    mutate(leadingEdge = map(leadingEdge, ~ .x[1:min(30, length(.x))])) %>%  # Limit to first 20 genes
    unnest(leadingEdge)
  
  significant_pathways <- filtered_pathways %>%
    distinct(pathway, NES) %>%
    arrange(desc(NES)) %>%
    pull(pathway)
  
  for (pathway in significant_pathways) {
    # Get genes for the current pathway
    pathway_genes <- filtered_pathways %>% 
      filter(pathway == !!pathway) %>% 
      pull(leadingEdge) %>% 
      unique()
    
    # Filter the matrix X for these genes
    pathway_X <- X[pathway_genes, c(sample1, sample2), drop = FALSE]
    
    # Skip pathway if no genes are found in X
    if (nrow(pathway_X) == 0) {
      next
    }
    
    pathway_info <- filtered_pathways %>%
      filter(pathway == !!pathway) %>%
      select(padj, NES) %>%
      distinct()
    
    padj_value <- signif(pathway_info$padj, 3)
    nes_value <- signif(pathway_info$NES, 3)
    
    # Define a color palette for the heatmap
    color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
    
    # Define breaks for the heatmap to limit range to [0, 5]
    breaks <- seq(0, 8, length.out = 101)
    
    # Create the heatmap for this pathway
    pheatmap(
      pathway_X,
      cluster_rows = FALSE,  # Cluster genes within the pathway
      cluster_cols = TRUE,  # Cluster columns (samples)
      color = color_palette,  # Use custom color palette
      # breaks = breaks,  # Set range to [0, 5]
      fontsize_row = 8,  # Adjust row font size
      fontsize_col = 8,  # Adjust column font size
      main = paste0("Pathway: ", pathway, "\n", "padj: ", padj_value, ", NES: ", nes_value),  # Add title
      border_color = NA  # Remove cell borders for a cleaner look
    )
    
  }
  dev.off()  # Close the PDF device
}


GSEA_diff <- function(expr_matrix, sample1, sample2, gene_list, GSEA_dir, n_permutations=10000, n=5) {
  log2_diff <- log2((expr_matrix[, sample2] + 1) / (expr_matrix[, sample1] + 1))
  
  ranked_genes <- sort(log2_diff, decreasing = TRUE)
  
  # gsea_results <- fgsea(pathways = gene_list,
  #                       stats = ranked_genes,
  #                       nperm = n_permutations)
  gsea_results <- fgseaMultilevel(pathways = gene_list,stats = ranked_genes)
  gsea_results$Log10padj <- -log10(gsea_results$padj)

  top_up <- gsea_results %>% arrange(desc(NES)) %>% slice_head(n = n)
  top_down <- gsea_results %>% arrange(NES) %>% slice_head(n = n)
  gsea_top <- bind_rows(top_up, top_down)
  gsea_top <- gsea_top %>% filter(padj < 0.05)
  gsea_sig <- gsea_results %>% filter(padj < 0.05)
  
  sample1N = sample1
  if (sample1N == "0") {
    sample1N = "R"
  }
  horizontal_top_plot <- ggplot(gsea_top, aes(x = Log10padj, y = reorder(pathway, NES), fill = NES)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = -log10(0.05), color = "green", linetype = "dashed", size = 1) +
    scale_fill_gradientn(colors = c("blue", "red"), name = "NES", oob = scales::squish, limits = c(-2, 2)) +
    labs(x = "-log10(p-adj)", y = "Pathway", title = paste("Clone ", sample2, " vs Clone ", sample1N)) +
    theme(axis.text.y = element_text(size = 24), 
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20),           # X-axis title size
          axis.title.y = element_text(size = 20),           # Y-axis title size
          plot.title = element_text(size = 26) # Title size and style
    )
  
  # Save the horizontal barplot for top pathways to a file
  output_file_top <- paste0(GSEA_dir,"/clone", sample2, "_", sample1, "_GSEA_diff", ".png") # Replace with desired file name
  ggsave(output_file_top, plot = horizontal_top_plot, width = 13, height = 12)
  
  
  horizontal_sig_plot <- ggplot(gsea_sig, aes(x = Log10padj, y = reorder(pathway, NES), fill = NES)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = -log10(0.05), color = "green", linetype = "dashed", size = 1) +
    scale_fill_gradientn(colors = c("blue", "red"), name = "NES", oob = scales::squish, limits = c(-2, 2)) +
    labs(x = "-log10(p-adj)", y = "Pathway", title = paste("Clone ", sample2, " vs Clone ", sample1N)) +
    theme(axis.text.y = element_text(size = 24), 
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20),           # X-axis title size
          axis.title.y = element_text(size = 20),           # Y-axis title size
          plot.title = element_text(size = 26) # Title size and style
    )
  
  # Save the horizontal barplot for top pathways to a file
  output_file_sig <- paste0(GSEA_dir,"/clone", sample2, "_", sample1, "_GSEA_sig", ".png") # Replace with desired file name
  ggsave(output_file_sig, plot = horizontal_sig_plot, width = 13, height = 12)
  
  gsea_df <- as.data.frame(gsea_results)
  gsea_df$leadingEdge <- sapply(gsea_df$leadingEdge, function(x) paste(x, collapse = ";"))
  write.csv(as.data.frame(gsea_df), file = paste0(GSEA_dir,"/clone", sample2, "_", sample1, "_GSEA_diff.csv"))
  # volcano_plot <- volcanoplot(gsea_results, paste0("Volcano plot for Top Differential Pathways between Clone ", sample2, " and Clone ", sample1), n=n)
  
  return(gsea_results)
}

read_tree <- function(treeFile) {
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
  
  return(edge_list)
}

volcanoplot <- function(data_full, title, n=5) {
  data_full <- data_full %>%
    mutate(Significance = ifelse(padj < 0.05, 
                                 ifelse(NES < 0, "Negative", "Positive"), "Not Significant"))
  
  labels <- data_full %>%
    arrange(desc(NES)) %>%
    slice_head(n = n) %>%
    bind_rows(
      data_full %>%
        arrange(NES) %>%
        slice_head(n = n)
    )
  
  volcano_plot <- ggplot(data_full, aes(x = NES, y = Log10padj, color = Significance)) +
    geom_point(size = 3, alpha = 0.8) + # Add points with size and transparency
    scale_color_manual(values = c("Negative" = "blue", "Positive" = "red", "Not Significant" = "gray")) + # Custom colors
    # geom_text(data = labels, aes(label = Pathway), vjust = -1, size = 3.5) + # Add labels
    geom_text_repel(data = labels, aes(label = pathway), size = 3.5, max.overlaps = Inf) + # Non-overlapping labels
    labs(
      title = title,
      x = "Enrichment Score",
      y = "-log10(P-value)",
      color = "Significant (p < 0.05)"
    ) +
    # xlim(-1, 1) + # Set x-axis limits
    # ylim(0, 3.2) +  # Set y-axis limits
    theme_minimal() +
    theme(legend.position = "top")
  
  return(volcano_plot)
}

read_GSEA_file <- function(GSEA_file) {
  lines <- readLines(GSEA_file) 
  
  gene_list <- list()
  
  for (line in lines) {
    elements <- strsplit(line, "\t")[[1]]
    gene_set_name <- elements[1]
    genes <- elements[-c(1,2)]
    gene_list[[gene_set_name]] <- genes
  }
  
  return(gene_list)
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
