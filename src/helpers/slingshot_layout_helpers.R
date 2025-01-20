# function wrappers to process slingshot outputs



################################################################################
################################################################################
pseudotime_branching <- function(pto,
                                 layout_spec, # either a function to layout the graph, or a colnames from sling_cell_data output
                                 layout_type = c(
                                   "default",
                                   "layout_only",
                                   "layout_dominant",
                                   "pseudotime_dominant",
                                   "pseudotime_only"
                                 ),
                                 offset = 10, 
                                 ...) {
  library(slingshot)
  library(RColorBrewer)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  
  
  # match options
  layout_type <- match.arg(layout_type)
  
  
  # retrieve cluster labels
  pseudotime <- sling_cell_data(pto)
  
  # choose cluster calls
  pseudotime <- pseudotime %>%
    mutate(group = as.character(cluster))
  
  # retrieve MST
  mst_graph <- slingMST(pto)
  
  root_node <- slingParams(pto)[["start.clus"]]
  
  # calculate average pseudotime data for segments
  pseudotime_means <- pseudotime %>%
    group_by(group) %>%
    summarize(pseudotime = mean(pseudotime))
  
  
  
  cat("root node", root_node, sep = "\n")
  
  
  
  # build nodes positions from the MST graph
  if(is.function(layout_spec)) {
    
    if(is.null(formals(layout_spec)[["root"]])) {
      
      graph_layout <- layout_spec(mst_graph, ...) %>%
        as.data.frame()
    } else {
      graph_layout <- layout_spec(mst_graph, root = root_node, ...) %>%
        as.data.frame()
    }
    
    colnames(graph_layout) <- c("x", "y")
    
  } else if(length(layout_spec) == 1) {
    
    graph_layout <- pseudotime %>%
      dplyr::select(x = layout_spec, id) %>%
      tibble::column_to_rownames("id")
    
  }
  
  
  
  
  
  nodes_position <- graph_layout %>%
    tibble::rowid_to_column("id") %>%
    mutate(group = vertex_attr(mst_graph)[["name"]]) %>%
    #mutate(group = unclass(as.factor((group)))) %>%
    left_join(pseudotime_means)
  
  
  
  
  # layout option to space the cluster on the branches
  
  if (layout_type == "layout_dominant") {
    nodes_position <- nodes_position %>%
      arrange(x, pseudotime)
  } else if (layout_type == "pseudotime_dominant") {
    nodes_position <- nodes_position %>%
      arrange(pseudotime, x)
  } else if (layout_type == "pseudotime_only") {
    nodes_position <- nodes_position %>%
      arrange(pseudotime)
  } else if (layout_type == "layout_only") {
    nodes_position <- nodes_position %>%
      arrange(x)
  }
  
  if (layout_type != "default") {
    nodes_position <- nodes_position %>%
      tibble::rowid_to_column("xx") %>%
      dplyr::mutate(x = x * offset + xx) %>%
      arrange(id)
  }
  
  
  
  
  # Assign cells to the node positions
  pseudotime <- pseudotime %>%
    left_join(nodes_position %>%
                dplyr::select(-c(pseudotime, id))) # add spacings
  
  
  plot <- ggraph(mst_graph, layout = nodes_position %>% dplyr::select(x = pseudotime, y = x), circular = FALSE) + 
          geom_edge_diagonal(flipped = T) +
          geom_node_label(aes(label = name)) +
          theme_bw() 
  
  
  
  plot +
    geom_violin(
      data = pseudotime,
      aes(
        x = pseudotime,
        y = x,
        fill = as.factor(cluster),
        #fill = after_stat(x),
        #height = after_stat(density)
      ),
      #stat = "density",
      position = "identity"
      #scale = 0.8,
      #trim = TRUE,
      
    ) +
    geom_label(data = nodes_position, aes(pseudotime, x, label = group))
  #geom_point(data = pseudotime_means, aes(x = pseudotime, y = x, color = group), shape = 1)
}

################################################################################
################################################################################
sling_cell_data <- function(pto) {
  
  library(slingshot)
  library(tidyverse)
  
  curve_probs <- slingCurveWeights(pto, as.probs = TRUE)
  
  #curve_choice <- curve_probs == rowMaxs(curve_probs)
  
  
  scaled_PT <- slingPseudotime(pto, na = F) * curve_probs
  
  
  scaled_PT[is.na(scaled_PT)]<-0
  colnames(scaled_PT) <- paste0("PT", 1:ncol(scaled_PT))
  
  pca <- prcomp(scaled_PT)[["x"]]
  
  avgP <- rowSums(scaled_PT)
  
  
  sling_data <- cbind(scaled_PT, pca) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    mutate(pseudotime = avgP) %>% # add pseudotime
    mutate(lineage = slingBranchID(pto))
  
  
  
  # retrieve cluster labels
  slingClusterLabels(pto) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "cluster", values_to = "keep") %>%
    dplyr::filter(keep == 1) %>%
    dplyr::select(-keep) %>%
    right_join(sling_data)
  
  
  
}


################################################################################
################################################################################
# does not work super well in putting the curves in the right place so far
pseudotime_lineages <- function(pto, smoothing = 0.5) {
  
  # retrieve cluster labels
  pseudotime <- sling_cell_data(pto)
  
  # if we want curve connector we'll fit the cells onto the curve + jittering
  # so fit the curves to the acyclic layout and plot the cells, and done
  
    edges <-  embedCurves(
      pto,
      dplyr::select(pseudotime, x = PC1, y = PC2) %>% as.matrix(),
      shrink = 0,
      #approx_points = 100,
      #smoother = "loess",
      #span = smoothing
    ) %>%
      slingCurves(as.df = T)
    
    print(edges)
    
    
    #return(edges)
    #print(edges)
    
    
    plot <- pseudotime %>%
      ggplot() +
      geom_line(data = edges, aes(x, y, group = as.factor(Lineage))) +
      geom_point(aes(PC1, PC2, color = as.factor(cluster)),
                 position = position_jitter(height = 0.5, width = 0.5, seed = 134))
    
    return(plot)
  
  
  
}

