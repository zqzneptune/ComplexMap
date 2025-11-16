utils::globalVariables(c("betweenness", "degree"))
#' Compute Network Topology and Layout Coordinates
#'
#' @description
#' Calculates the 2D layout coordinates and key centrality metrics for a complex
#' network. This function serves as the final step in preparing network data
#' for visualization.
#'
#' @details
#' This function takes a node attribute table and an edge list (network) and
#' performs the following steps:
#' 
#' 1.  Constructs an `igraph` graph object from the provided data.
#' 
#' 2.  Computes a force-directed layout using the Fruchterman-Reingold algorithm
#'     via `ggraph::create_layout`. Edge weights are used to influence the
#'     layout, pulling strongly connected nodes closer together.
#'     
#' 3.  Calculates node centrality metrics:
#' 
#'     - **Betweenness Centrality:** Measures how often a node lies on the
#'       shortest path between other nodes (normalized).
#'       
#'     - **Degree Centrality:** The number of edges connected to a node.
#'     
#' 4.  Merges the layout coordinates and centrality scores back into the original
#'     node attribute table.
#'
#' @param nodeAttributes A `tibble` or `data.frame` containing attributes for
#'   each node (complex). Must contain a column with node identifiers that
#'   matches the source/target columns in the `network` data.
#' @param network A `tibble` or `data.frame` representing the network edges.
#'   Must contain columns for source, target, and edge `weight`.
#' @param verbose A logical value indicating whether to print progress messages.
#'   Defaults to `TRUE`.
#'
#' @return
#' A `tibble` containing all original columns from `nodeAttributes` plus four
#' new columns: `x`, `y` (layout coordinates), `betweenness`, and `degree`.
#' The table is arranged in descending order of betweenness and degree.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @seealso
#' `generateNodeAttributes()`, `buildComplexNetwork()`
#'
#' @examples
#' # --- Sample Data ---
#' # 1. Node attributes
#' nodes <- tibble::tibble(
#'   complexId = c("Cpx1", "Cpx2", "Cpx3"),
#'   proteinCount = c(10, 8, 12)
#' )
#'
#' # 2. Network edges
#' net <- tibble::tibble(
#'   source = c("Cpx1", "Cpx2"),
#'   target = c("Cpx2", "Cpx3"),
#'   weight = c(0.8, 0.6)
#' )
#'
#' # --- Compute Topology ---
#' masterLayout <- computeMapTopology(nodes, net)
#' print(masterLayout)
#'
#' @export
#'
computeMapTopology <- function(nodeAttributes, network, verbose=TRUE) {
  if (verbose) {
    message("Computing map topology (layout and centrality)...")
  }
  
  # Handle cases with empty input data gracefully
  if (nrow(nodeAttributes) == 0 || nrow(network) == 0) {
    warning("Cannot compute topology with empty node or edge lists.")
    # Return a dataframe with the right columns but no data
    return(
      nodeAttributes %>%
        dplyr::mutate(
          x = NA_real_, y = NA_real_,
          betweenness = NA_real_, degree = NA_integer_
        )
    )
  }
  
  # Create the graph object, ensuring the node ID column is correctly identified
  # igraph uses the first column of the vertices df as the name by default.
  # Assuming the first column of nodeAttributes holds the complex IDs.
  graphObj <- igraph::graph_from_data_frame(
    d=network, vertices=nodeAttributes, directed=FALSE
  )
  
  # 1. Compute the force-directed layout
  layoutData <- ggraph::create_layout(
    graphObj, layout='fr', weights=igraph::E(graphObj)$weight
  )
  
  # 2. Compute centrality metrics
  betweennessScores <- igraph::betweenness(graphObj, normalized=TRUE)
  degreeScores <- igraph::degree(graphObj)
  
  # 3. Combine everything into the final master dataframe
  # The first column of nodeAttributes must be the complex_id for this to work
  idColName <- names(nodeAttributes)[1]
  
  masterLayoutDf <- nodeAttributes %>%
    dplyr::mutate(
      x = layoutData$x,
      y = layoutData$y,
      betweenness = betweennessScores[!!rlang::sym(idColName)],
      degree = degreeScores[!!rlang::sym(idColName)]
    ) %>%
    dplyr::arrange(dplyr::desc(betweenness), dplyr::desc(degree))
  
  if (verbose) {
    message("Topology computation complete.")
  }
  
  return(masterLayoutDf)
}