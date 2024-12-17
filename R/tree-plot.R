#' Plot ensemble tree
#' 
#' @export
#' @param trees list of tibbles of edges, each with columns edge, parent, child
plotEnsembleTree <- function(trees, palette=viridis::viridis) {
  am_chain <- lapply(trees, edgesToAmLong)
  post_am <- getPosteriorAmLong(am_chain)
  plotPosteriorAmLong(post_am, colorScheme(trees[[1]], palette))
}

plotPosteriorAmLong <- function(post_am, v_color, filter1 = TRUE, filter1.threshold = 0.1,
                                filter2 = TRUE, filter2.threshold = 0.1) {
  # filter1 filters columns (am wide format) for edges with posterior prob > (max(column) - filter1.threshold)
  admat <- prepPostAmForGraphing(post_am)
  
  # filter edges of low freq
  admat <- filterAdmat(admat, filter1 = filter1, filter1.threshold = filter1.threshold,
                       filter2 = filter2, filter2.threshold = filter2.threshold)
  
  ig <- igraph::graph_from_adjacency_matrix(admat, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE) 
  
  igraph::E(ig)$lty <- ifelse(igraph::E(ig)$weight < 0.25, 2, 1)
  
  # make edge black if only 1 edge to vertex
  e <- igraph::ends(ig, igraph::E(ig))
  numTo <- table(e[,2])
  edgeColors <- sapply(e[,2], function(x) ifelse(x %in% names(which(numTo==1)), "black", "darkgrey"))
  igraph::E(ig)$color <- edgeColors
  
  igraph::V(ig)$label.cex <- 1.5
  
  igraph::V(ig)$color <- as.list(v_color %>% arrange(match(v_sorted, names(V(ig)))) %>% select(colors))$colors
  
  par(mar=c(0,0,0,0)+.1)
  
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      vertex.label.family = "Helvetica", vertex.size=20,
                      edge.arrow.size = 0.5, edge.arrow.width = 2,
                      edge.width = igraph::E(ig)$weight*3)
}

#' Plot single tree 
#' 
#' @export
#' @param edges tibble of edges with columns edge, parent, child
plotTree <- function(edges, filtered_table=NULL, palette=viridis::viridis) {
  # plotGraph(edgesToAmLong(edges), colorScheme(edges, palette), filtered_table)
  
  v_color <- colorScheme(edges, palette)
  edge_list <- as.matrix(edges[, c("parent", "child")])
  ig <- graph_from_edgelist(edge_list, directed = TRUE)
  # E(ig)$name <- edges$edge
  ###########
  
  # ig <- igraph::graph_from_adjacency_matrix(am, mode = "directed", weighted = TRUE,
  #                                           diag = FALSE, add.row = TRUE)
  V(ig)$color <- as.list(v_color %>% arrange(match(v_sorted, names(V(ig)))) %>% select(colors))$colors
  
  V(ig)$label.color <- ifelse(calculate_brightness(V(ig)$color) < 100, "white", "black")
  
  par(mar=c(1,1,1,1))
  
  ##################
  ig <- set_edge_attr(ig, "Mut_ID", value = "")
  
  if (!is.null(filtered_table)) {
    
    for (i in seq_len(nrow(filtered_table))) {
      cluster_name <- filtered_table$Cluster[i]
      new_mut_id <- filtered_table$Mut_ID[i]
      
      # Find the node index for the matching cluster
      node_index <- which(V(ig)$name == cluster_name)
      
      if (length(node_index) > 0) {
        # Find the edge leading to this node
        in_edge <- incident(ig, node_index, mode = "in") # Incoming edge to the node
        
        if (length(in_edge) > 0) {
          # Append the new_Mut_ID to the edge attribute
          current_label <- edge_attr(ig, "Mut_ID", index = in_edge)
          updated_label <- ifelse(
            current_label == "",
            new_mut_id,
            paste(current_label, new_mut_id, sep = "\n ")
          )
          ig <- set_edge_attr(ig, "Mut_ID", index = in_edge, value = updated_label)
        }
      }
    }
  }
  
  edge_lengths <- sapply(E(ig)$Mut_ID, function(label) {
    if (label == "") {
      1  # Default length if no label
    } else {
      str_count(label, "\n") + 1  # Count lines (number of \n + 1)
    }
  })
  
  min_length <- 3  # Adjust this value as needed
  adjusted_edge_lengths <- pmax(edge_lengths, min_length)
  
  layout <- layout_as_tree(ig)
  
  max_cluster <- filtered_table %>%
    count(Cluster) %>%                  # Count occurrences of each unique Cluster
    arrange(desc(n)) %>%                # Arrange by count in descending order
    slice(1) 
  
  if (ecount(ig) > 3 & max_cluster$n > 4) {
    # Scale the y-coordinates of the layout based on edge lengths
    for (i in seq_len(ecount(ig))) {
      parent <- ends(ig, i)[1]  # Get the parent node
      child <- ends(ig, i)[2]   # Get the child node
      
      # Scale child y-coordinate based on edge length
      layout[which(V(ig)$name == child), 2] <- 
        layout[which(V(ig)$name == parent), 2] - adjusted_edge_lengths[i] * 0.1
    }
  }
 
  
  ########
  V(ig)$name[V(ig)$name == "root"] <- "R"
  
  igraph::plot.igraph(ig, layout = layout, edge.label = edge_attr(ig, "Mut_ID"),
                      vertex.size=26, vertex.frame.color = "#000000", vertex.label.cex = 1.2, edge.label.cex=0.5,
                      vertex.label.family = "Helvetica", vertex.label.color = V(ig)$label.color, 
                      edge.arrow.size = 0.4, edge.arrow.width = 1, edge.color = "grey", edge.label.color = "red",
                      margin=0.2)
  
  # igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig), edge.label = edge_attr(ig, "Mut_ID"),
  #                     vertex.size=32, vertex.frame.color = "#000000", vertex.label.cex = 1.6, edge.label.cex=0.75,
  #                     vertex.label.family = "Helvetica", vertex.label.color = V(ig)$label.color, 
  #                     edge.arrow.size = 0.5, edge.arrow.width = 2, edge.color = "black", edge.label.color = "red",
  #                     margin=0.2)
  
}

#' generate colors for each vertice
#' @export
colorScheme <- function(edges, palette=viridis::viridis) {
  v_sorted = sort(unique(c(edges$parent, edges$child)))
  v_sorted = c(sort(as.integer(v_sorted[!v_sorted=='root'])), "root")
  # root_idx <- which(v_sorted=="root")
  colors <- c(palette(length(v_sorted)-1), "white")
  v_color <- tibble(v_sorted, colors)
  return(v_color)
}

calculate_brightness <- function(hex_color) {
  rgb <- col2rgb(hex_color)
  # Brightness formula: 0.299*R + 0.587*G + 0.114*B
  brightness <- 0.299 * rgb[1,] + 0.587 * rgb[2,] + 0.114 * rgb[3,]
  return(brightness)
}

plotGraph <- function(am.long, v_color, filtered_table=NULL){
  # make sure am.long is sorted by parent and child
  am.long <- mutate(am.long, child = as.numeric(am.long$child)) %>%
    arrange(parent, child)
  am.long <- mutate(am.long, child = as.character(am.long$child))
  
  # change to wide format and plot
  am <- toWide(am.long)
  rownames(am) <- c("root", colnames(am))
  am <- cbind(root=0, am) ## add column for root
  colnames(am) <- rownames(am)
  
  am[is.na(am)] <- 0
  
  ig <- igraph::graph_from_adjacency_matrix(am, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE)
  V(ig)$color <- as.list(v_color %>% arrange(match(v_sorted, names(V(ig)))) %>% select(colors))$colors
  
  V(ig)$label.color <- ifelse(calculate_brightness(V(ig)$color) < 100, "white", "black")
  
  par(mar=c(3,3,3,3))
  
  ##################
  ig <- set_edge_attr(ig, "Mut_ID", value = "")
  
  if (!is.null(filtered_table)) {
    
    for (i in seq_len(nrow(filtered_table))) {
      cluster_name <- filtered_table$Cluster[i]
      new_mut_id <- filtered_table$Mut_ID[i]
      
      # Find the node index for the matching cluster
      node_index <- which(V(ig)$name == cluster_name)
      
      if (length(node_index) > 0) {
        # Find the edge leading to this node
        in_edge <- incident(ig, node_index, mode = "in") # Incoming edge to the node
        
        if (length(in_edge) > 0) {
          # Append the new_Mut_ID to the edge attribute
          current_label <- edge_attr(ig, "Mut_ID", index = in_edge)
          updated_label <- ifelse(
            current_label == "",
            new_mut_id,
            paste(current_label, new_mut_id, sep = "\n ")
          )
          ig <- set_edge_attr(ig, "Mut_ID", index = in_edge, value = updated_label)
        }
      }
    }
  }
  
  ########
  # V(ig)$name[V(ig)$name == "root"] <- "R"
  
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig), edge.label = edge_attr(ig, "Mut_ID"),
                      vertex.size=32, vertex.frame.color = "#000000", vertex.label.cex = 1.6, edge.label.cex=0.75,
                      vertex.label.family = "Helvetica", vertex.label.color = V(ig)$label.color,
                      edge.arrow.size = 0.5, edge.arrow.width = 2, edge.color = "black", edge.label.color = "red",
                      margin=0.3)
  
  

}