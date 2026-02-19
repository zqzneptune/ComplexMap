utils::globalVariables(c("x", "y", "weight", "source_complex_id", "target_complex_id",
                         "compSim", "funcSim", "complexId", "proteinCount",
                         "primaryFunctionalDomain", "colorHex", "betweenness",
                         "degree", "topEnrichedFunctions", "proteins"))

# ============================================================================
# CYTOSCAPE.JS EXPORT LAYER
# ============================================================================

#' Scale Layout Coordinates for Cytoscape.js
#'
#' @description
#' Mean-centers and rescales raw layout coordinates to a fixed pixel range
#' suitable for Cytoscape.js.
#'
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates.
#' @param range Numeric vector of length 2 giving target range. Default `c(-500, 500)`.
#'
#' @return A `list` with elements `x` and `y`, both rescaled.
#'
#' @export
scaleLayoutForCytoscape <- function(x, y, range = c(-500, 500)) {
  scale_vec <- function(v, range) {
    mn  <- min(v, na.rm = TRUE)
    mx  <- max(v, na.rm = TRUE)
    if (mn == mx) return(rep((range[1] + range[2]) / 2, length(v)))
    range[1] + (v - mn) / (mx - mn) * (range[2] - range[1])
  }
  list(x = scale_vec(x, range), y = scale_vec(y, range))
}

#' Prune Network Edges Before Cytoscape Export
#'
#' @description
#' Reduces the number of edges in a `ComplexMap` object to prevent browser
#' overload during interactive exploration.
#'
#' @param cmap A `ComplexMap` object.
#' @param method Pruning method: `"top_k"` (keep top-k edges per node) or
#'   `"quantile"` (keep edges above a weight quantile).
#' @param k Integer. Number of top edges to retain per node when `method = "top_k"`.
#' @param weight_quantile Numeric (0-1). Quantile threshold when `method = "quantile"`.
#'   Defaults to `0.75`.
#' @param verbose Logical.
#'
#' @return A modified `ComplexMap` object with pruned edges.
#'
#' @export
pruneNetwork <- function(cmap, method = "top_k", k = 5,
                         weight_quantile = 0.75, verbose = TRUE) {
  if (!inherits(cmap, "ComplexMap")) {
    stop("'cmap' must be a ComplexMap object.", call. = FALSE)
  }

  edges <- getEdgeTable(cmap)
  original_n <- nrow(edges)

  if (original_n == 0) {
    if (verbose) message("No edges to prune.")
    return(cmap)
  }

  if (method == "quantile") {
    threshold <- stats::quantile(edges$weight, probs = weight_quantile, na.rm = TRUE)
    edges <- edges[edges$weight >= threshold, ]

  } else if (method == "top_k") {
    # For each node, keep its top-k neighbours by weight
    # Build adjacency list, keep top k per node, take union of all kept edges
    all_nodes <- unique(c(edges$source_complex_id, edges$target_complex_id))

    keep_idx <- lapply(all_nodes, function(node) {
      node_edges <- which(edges$source_complex_id == node | edges$target_complex_id == node)
      if (length(node_edges) <= k) return(node_edges)
      node_edges[order(edges$weight[node_edges], decreasing = TRUE)[seq_len(k)]]
    })

    keep_set <- sort(unique(unlist(keep_idx)))
    edges <- edges[keep_set, ]

  } else {
    stop("'method' must be 'top_k' or 'quantile'.", call. = FALSE)
  }

  if (verbose) {
    message(sprintf("pruneNetwork: %d -> %d edges (method='%s')", original_n, nrow(edges), method))
  }

  cmap$edges <- edges
  return(cmap)
}

#' Export a ComplexMap Object to Cytoscape.js JSON
#'
#' @description
#' Converts a `ComplexMap` object into a Cytoscape.js-compatible JSON string,
#' with scaled layout positions. This JSON is consumed directly by the
#' Shiny/Cytoscape.js interactive explorer.
#'
#' @param cmap A `ComplexMap` object (with layout coordinates in the node table).
#' @param prune Logical. Whether to prune the network before export. Default `TRUE`.
#' @param prune_method Pruning method passed to `pruneNetwork()`. Default `"top_k"`.
#' @param k Integer. Top-k edges per node if `prune_method = "top_k"`. Default `5`.
#' @param weight_quantile Numeric. Quantile threshold if `prune_method = "quantile"`. Default `0.75`.
#' @param coord_range Numeric length-2 vector for coordinate scaling. Default `c(-500, 500)`.
#' @param verbose Logical.
#'
#' @return A JSON string suitable for `window.complexmapElements` in JavaScript.
#'
#' @export
toCytoscapeJSON <- function(cmap,
                            prune          = TRUE,
                            prune_method   = "top_k",
                            k              = 5,
                            weight_quantile = 0.75,
                            coord_range    = c(-500, 500),
                            verbose        = TRUE) {

  if (!inherits(cmap, "ComplexMap")) {
    stop("'cmap' must be a ComplexMap object.", call. = FALSE)
  }

  # 1. Optionally prune edges
  if (prune) {
    cmap <- pruneNetwork(cmap, method = prune_method, k = k,
                         weight_quantile = weight_quantile, verbose = verbose)
  }

  nodes <- getNodeTable(cmap)
  edges <- getEdgeTable(cmap)

  # 2. Check required layout columns exist
  if (!all(c("x", "y") %in% names(nodes))) {
    stop("Node table lacks layout coordinates. Run computeMapTopology() first.", call. = FALSE)
  }

  # 3. Scale layout
  scaled <- scaleLayoutForCytoscape(nodes$x, nodes$y, range = coord_range)

  # 4. Build node elements
  node_elements <- lapply(seq_len(nrow(nodes)), function(i) {
    nd <- nodes[i, ]

    # Gather optional fields safely
    spec_score  <- if ("score"       %in% names(nd)) nd$score       else NA_real_
    top_funcs   <- if ("topEnrichedFunctions" %in% names(nd)) nd$topEnrichedFunctions else NA_character_

    list(
      data = list(
        id                      = nd$complexId,
        label                   = nd$complexId,
        proteinCount            = nd$proteinCount,
        primaryFunctionalDomain = nd$primaryFunctionalDomain,
        specificityScore        = if (is.na(spec_score)) 0 else round(spec_score, 4),
        degree                  = if ("degree"      %in% names(nd)) nd$degree                   else 0L,
        betweenness             = if ("betweenness" %in% names(nd)) round(nd$betweenness, 6)    else 0,
        colorHex                = nd$colorHex,
        proteins                = nd$proteins,
        topEnrichedFunctions    = if (is.na(top_funcs)) "" else top_funcs
      ),
      position = list(x = round(scaled$x[i], 2), y = round(scaled$y[i], 2))
    )
  })

  # 5. Build edge elements
  edge_elements <- lapply(seq_len(nrow(edges)), function(i) {
    ed <- edges[i, ]
    comp_sim <- if ("compSim" %in% names(ed)) round(ed$compSim, 4) else NA_real_
    func_sim <- if ("funcSim" %in% names(ed)) round(ed$funcSim, 4) else NA_real_

    list(
      data = list(
        id      = paste0("e", i),
        source  = ed$source_complex_id,
        target  = ed$target_complex_id,
        weight  = round(ed$weight, 4),
        compSim = comp_sim,
        funcSim = func_sim
      )
    )
  })

  all_elements <- c(node_elements, edge_elements)

  json_str <- jsonlite::toJSON(all_elements, auto_unbox = TRUE, na = "null", digits = 6)

  if (verbose) {
    message(sprintf("toCytoscapeJSON: %d nodes, %d edges exported.",
                    length(node_elements), length(edge_elements)))
  }

  return(json_str)
}
