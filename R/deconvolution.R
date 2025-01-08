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
  
  # edges <- list()
  # con <- file(treeFile, "r")
  # on.exit(close(con))
  # readLines(con, n = 1)
  # 
  # while(TRUE) {
  #   line <- readLines(con, n = 1)
  #   if (length(line) == 0) break  # Exit loop if no more lines
  #   line <- strsplit(line, ",")[[1]]
  #   
  #   # Add edges based on conditions
  #   if (line[2] == "root") {
  #     edges <- append(edges, list(c(0, as.integer(line[3]))))
  #   } else if (line[3] == "root") {
  #     edges <- append(edges, list(c(as.integer(line[2]), 0)))
  #   } else {
  #     edges <- append(edges, list(c(as.integer(line[2]), as.integer(line[3]))))
  #   }
  # }
  # 
  # # Create a graph and add edges
  # edge_list <- as.matrix(do.call(rbind, edges))  # Convert to matrix if not already
  # edge_list <- apply(edge_list, 2, as.integer) 
  
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
  
  # if (GSEA) {
  #   GSEA_analysis(X_optimal, outputDir, GSEA_file=GSEA_file)
  # }
  return(X_optimal)
  
}


#' @export
#' @import GSVA pheatmap limma ggplot2 ggrepel fgsea GSEABase
runGSEA <- function(X_optimal,
                    outputDir,
                    treeFile,
                    GSEA_file=NULL,
                    top_K=5,
                    n_permutations=1000) {
  
  GSEA_dir = paste0(outputDir, "/GSEA")
  suppressWarnings(dir.create(GSEA_dir))
  ssGSEA_dir = paste0(GSEA_dir, "/ssGSEA")
  suppressWarnings(dir.create(ssGSEA_dir))
  
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
  
  ############# ssGSEA2 #########
  # gctFile = paste0(GSEA_dir, "/clone_expression.gct")
  # writeGCT(X, gctFile)
  #   
  # res = run_ssGSEA2(gctFile,
  #                   output.prefix = "patient",
  #                   gene.set.databases = GSEA_file,
  #                   output.directory = ssGSEA_dir,
  #                   output.score.type = "NES", 
  #                   global.fdr = TRUE,
  #                   log.file = paste0(GSEA_dir, "/ssGSEA/run.log"))
  # 
  # 
  # ssGSEA_res <- read.delim(paste0(ssGSEA_dir, "/patient-combined.gct"), skip=2, header = TRUE)
  # 
  # sample_name = "1"
  # 
  # pvalue_col = paste0("pvalue.", sample_name)
  # fdr_col = paste0("fdr.pvalue.", sample_name)
  # nes_col = paste0("X", sample_name)
  # 
  # temp <- ssGSEA_res %>% 
  #   select(id, !!sym(pvalue_col), !!sym(fdr_col), !!sym(nes_col)) %>% 
  #   rename (pvalue = !!sym(pvalue_col), fdr = !!sym(fdr_col), NES =!!sym(nes_col))
  
  ############################
  
  # # run ssGSEA for each clone
  # gsvaPar <- ssgseaParam(X, gene_list)
  # gsva.es <- gsva(gsvaPar)
  # 
  # # permutation test to generate p values
  # perm_results <- matrix(0, nrow = nrow(gsva.es), ncol = ncol(gsva.es))
  # # n_permutations <- 10000
  # 
  # for (i in 1:n_permutations) {
  #   # Permute gene labels
  #   permuted_data <- X[sample(nrow(X)),,drop=FALSE]
  #   rownames(permuted_data) <- rownames(X)
  #   
  #   # Perform ssGSEA on permuted data
  #   permPar <- ssgseaParam(permuted_data, gene_list)
  #   perm_ssgsea <- suppressMessages(gsva(permPar))
  #   
  #   # Store results
  #   perm_results <- perm_results + ((gsva.es > 0 & gsva.es <= perm_ssgsea) | (gsva.es < 0 & gsva.es >= perm_ssgsea))
  # }
  # 
  # p_values <- perm_results / n_permutations
  # 
  # p_values_corrected <- matrix(p.adjust(as.vector(p_values), method = "fdr"),
  #                              nrow = nrow(p_values), 
  #                              ncol = ncol(p_values))
  # colnames(p_values_corrected) = colnames(p_values)
  # rownames(p_values_corrected) = rownames(p_values)
  
  # plotting for clone level GSEA
  
  ### ssGSEA using fgsea
  
  # samples <- colnames(X)
  # 
  # plots <- lapply(samples, function(sample_name) {
  #   
  #   sample_gsea <- ssGSEA(X, sample_name, gene_list, GSEA_dir, n_permutations=n_permutations, n=top_K)
  #   # 
  #   # data_full <- prepare_barplot_data(gsva.es, p_values_corrected, sample_name)
  #   # 
  #   # write.csv(data_full, file = paste0(GSEA_dir,"/clone", sample_name, "_ssGSEA.csv"))
  #   # 
  #   # data_top <- filter_top_pathways(data_full, n=top_K)
  #   # 
  #   # # Top pathways barplot
  #   # horizontal_top_plot <- ggplot(data_top, aes(x = Log10Pval, y = reorder(Pathway, Enrichment), fill = Enrichment)) +
  #   #   geom_bar(stat = "identity") +
  #   #   geom_vline(xintercept = -log10(0.05), color = "green", linetype = "dashed", size = 1) +
  #   #   scale_fill_gradientn(colors = c("blue", "red"), name = "NES", oob = scales::squish) +
  #   #   labs(x = "-log10(p-adj)", y = "Pathway", title = paste("Top Enrichment Pathways for Clone ", sample_name)) +
  #   #   theme(axis.text.y = element_text(size = 16), 
  #   #         axis.text.x = element_text(size = 14),
  #   #         axis.title.x = element_text(size = 18),           # X-axis title size
  #   #         axis.title.y = element_text(size = 18),           # Y-axis title size
  #   #         plot.title = element_text(size = 20) # Title size and style
  #   #   )
  #   # 
  #   # # Save the horizontal barplot for top pathways to a file
  #   # output_file_top <- paste0(GSEA_dir,"/clone", sample_name, "_ssGSEA_top", ".png") # Replace with desired file name
  #   # ggsave(output_file_top, plot = horizontal_top_plot, width = 13, height = 12)
  #   # 
  #   # volcano plot
  #   
  # })
  
  edge_list <- read_tree(treeFile)
  
  gsea_results_list <- list()
  
  if (is.null(ncol(edge_list))) {
    
    edge_list <- as.character(edge_list)
    
    sample1 = edge_list[1]
    sample2 = edge_list[2]
    
    gsea_results_list[[1]] <- GSEA_diff(X, sample1, sample2, gene_list, GSEA_dir, n_permutations, top_K)
    # gsea_results <- GSEA_diff(X, sample1, sample2, gene_list, n_permutations)
    
  } else {
    edge_list <- apply(edge_list, 2, as.character)
    
    for (i in 1:nrow(edge_list)) {
      
      sample1 = edge_list[i, 1]
      sample2 = edge_list[i, 2]
      
      gsea_results_list[[i]] <- GSEA_diff(X, sample1, sample2, gene_list, GSEA_dir, n_permutations, top_K)
      # gsea_results <- GSEA_diff(X, sample1, sample2, gene_list, n_permutations)
      
    }
  }
  
  # length(gsea_results_list)
  
}

writeGCT <- function(gene_counts, file) {
  
  gct_data <- as.data.frame(gene_counts)
  # gct_data$Description <- "NA"
  gct_data <- cbind(id = rownames(gct_data), gct_data)
  rownames(gct_data) <- NULL
  
  cat("#1.3\n", file = file)
  cat(nrow(gct_data), ncol(gct_data) - 1, 0, 0, "\n", file = file, append = TRUE, sep="\t")
  
  # Write data
  write.table(gct_data, file = file, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
  
  cat("GCT file written to:", file, "\n")
}

ssGSEA <- function(X, sample_name, gene_list, GSEA_dir, n_permutations=10000, n=5) {
  sample_ranked_genes <- sort(X[,sample_name], decreasing = T)
  
  # sample_gsea <- fgsea(pathways = gene_list, stats = sample_ranked_genes, nperm = n_permutations)
  sample_gsea <- fgseaMultilevel(pathways = gene_list, stats = sample_ranked_genes)
  
  
  sample_gsea$Log10padj <- -log10(sample_gsea$padj)
  
  # top_up <- sample_gsea %>% arrange(desc(NES)) %>% slice_head(n = n)
  # top_down <- sample_gsea %>% arrange(NES) %>% slice_head(n = n)
  # sample_gsea_top <- bind_rows(top_up, top_down)
  sample_gsea_top <- sample_gsea %>% filter(padj < 0.05)
  
  horizontal_top_plot <- ggplot(sample_gsea_top, aes(x = Log10padj, y = reorder(pathway, NES), fill = NES)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = -log10(0.05), color = "green", linetype = "dashed", linewidth = 1) +
    scale_fill_gradientn(colors = c("blue", "red"), name = "NES", oob = scales::squish) +
    labs(x = "-log10(p-adj)", y = "Pathway", title = paste0("Top Enrichment Pathways for Clone ", sample_name)) +
    theme(axis.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),           # X-axis title size
          axis.title.y = element_text(size = 18),           # Y-axis title size
          plot.title = element_text(size = 20) # Title size and style
    )
  
  # Save the horizontal barplot for top pathways to a file
  output_file_top <- paste0(GSEA_dir,"/clone", sample_name, "_ssGSEA_top", ".png") # Replace with desired file name
  ggsave(output_file_top, plot = horizontal_top_plot, width = 13, height = 12)
  
  volcano_plot <- volcanoplot(sample_gsea, paste0("Volcano plot for Clone ", sample_name), n=n)
  
  output_file_top <- paste0(GSEA_dir,"/clone", sample_name, "_ssGSEA_top_volcano", ".png") # Replace with desired file name
  ggsave(output_file_top, plot = volcano_plot, width = 13, height = 12)
  
  leadingEdge <- sample_gsea$leadingEdge
  sample_gsea$leadingEdge <- NULL
  write.csv(as.data.frame(sample_gsea), file = paste0(GSEA_dir,"/clone", sample_name, "_ssGSEA.csv"))
  sample_gsea$leadingEdge <- leadingEdge
  return(sample_gsea)
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
  # gsea_top <- gsea_results %>% filter(padj < 0.05)
  
  horizontal_top_plot <- ggplot(gsea_top, aes(x = Log10padj, y = reorder(pathway, NES), fill = NES)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = -log10(0.05), color = "green", linetype = "dashed", size = 1) +
    scale_fill_gradientn(colors = c("blue", "red"), name = "NES", oob = scales::squish) +
    labs(x = "-log10(p-adj)", y = "Pathway", title = paste("Top Differential Pathways between Clone ", sample2, " and Clone ", sample1)) +
    theme(axis.text.y = element_text(size = 16), 
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),           # X-axis title size
          axis.title.y = element_text(size = 18),           # Y-axis title size
          plot.title = element_text(size = 20) # Title size and style
    )
  
  # Save the horizontal barplot for top pathways to a file
  output_file_top <- paste0(GSEA_dir,"/clone", sample2, "_", sample1, "_GSEA_diff", ".png") # Replace with desired file name
  ggsave(output_file_top, plot = horizontal_top_plot, width = 13, height = 12)
  
  leadingEdge <- gsea_results$leadingEdge
  gsea_results$leadingEdge <- NULL
  write.csv(as.data.frame(gsea_results), file = paste0(GSEA_dir,"/clone", sample2, "_", sample1, "_GSEA_diff.csv"))
  gsea_results$leadingEdge <- leadingEdge
  
  volcano_plot <- volcanoplot(gsea_results, paste0("Volcano plot for Top Differential Pathways between Clone ", sample2, " and Clone ", sample1), n=n)

  output_file_top <- paste0(GSEA_dir,"/clone", sample2, "_", sample1, "_GSEA_diff_volcano", ".png") # Replace with desired file name
  ggsave(output_file_top, plot = volcano_plot, width = 13, height = 12)
  
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

prepare_barplot_data <- function(enrichment_matrix, pval_matrix, sample_name) {
  enrichment <- enrichment_matrix[, sample_name]
  pval <- pval_matrix[, sample_name]
  
  # Create a data frame
  df <- data.frame(
    Pathway = rownames(enrichment_matrix),
    Enrichment = enrichment,
    Log10Pval = -log10(pval)
  )
  
  df$Log10Pval[df$Log10Pval == Inf] <- 3
  
  # Rank pathways by enrichment score
  df <- df %>% arrange(desc(Enrichment)) %>% mutate(Rank = row_number())
  return(df)
}

filter_top_pathways <- function(df, n=5) {
  top_up <- df %>% arrange(desc(Enrichment)) %>% slice_head(n = n)
  top_down <- df %>% arrange(Enrichment) %>% slice_head(n = n)
  return(bind_rows(top_up, top_down))
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
