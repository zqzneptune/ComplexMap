#' Run the ComplexMap Interactive Explorer
#'
#' @description
#' Launches a Shiny application that renders a `ComplexMap` object as an
#' interactive Cytoscape.js network. Nodes are colored by functional domain,
#' sized by protein count, and fully filterable via sidebar controls.
#'
#' @param cmap A `ComplexMap` object.
#' @param prune Logical. Whether to prune edges before launching. Default `TRUE`.
#' @param prune_method Pruning method: `"top_k"` or `"quantile"`. Default `"top_k"`.
#' @param k Integer. Top-k edges per node when `prune_method = "top_k"`. Default `5`.
#' @param weight_quantile Numeric (0-1). Quantile threshold for `"quantile"` method. Default `0.75`.
#' @param port Optional integer port for the Shiny server.
#' @param launch.browser Logical. Open browser automatically. Default `TRUE`.
#' @param verbose Logical.
#'
#' @return Launches the Shiny app (invisibly returns `NULL`).
#'
#' @export
runComplexMapApp <- function(cmap,
                             prune          = TRUE,
                             prune_method   = "top_k",
                             k              = 5,
                             weight_quantile = 0.75,
                             port           = NULL,
                             launch.browser = TRUE,
                             verbose        = TRUE) {

  if (!inherits(cmap, "ComplexMap")) {
    stop("'cmap' must be a ComplexMap object.", call. = FALSE)
  }

  if (verbose) message("Preparing ComplexMap for interactive exploration...")

  # Convert to Cytoscape JSON
  elements_json <- toCytoscapeJSON(
    cmap,
    prune           = prune,
    prune_method    = prune_method,
    k               = k,
    weight_quantile = weight_quantile,
    verbose         = verbose
  )

  # Gather metadata for Shiny UI dropdowns
  nodes <- getNodeTable(cmap)
  edges <- getEdgeTable(cmap)

  domains <- sort(unique(nodes$primaryFunctionalDomain))

  # Pass data to Shiny via options (app-global)
  shiny::shinyOptions(
    complexmap_elements    = elements_json,
    complexmap_nodes       = nodes,
    complexmap_edges       = edges,
    complexmap_domains     = domains,
    complexmap_layout_info = cmap$layout_info
  )

  app_dir <- system.file("shiny", package = "ComplexMap")
  if (nchar(app_dir) == 0 || !dir.exists(app_dir)) {
    stop("Shiny app directory not found in ComplexMap package install.", call. = FALSE)
  }

  shiny::runApp(
    appDir         = app_dir,
    port           = port,
    launch.browser = launch.browser,
    quiet          = !verbose
  )

  invisible(NULL)
}

#' Explore a ComplexMap Object Interactively
#'
#' @description
#' A convenient alias for `runComplexMapApp()`. Launches the Cytoscape.js
#' interactive explorer.
#' @param ... Additional arguments passed to `runComplexMapApp`.
#'
#' @inheritParams runComplexMapApp
#'
#' @export
explore <- function(cmap, ...) {
  runComplexMapApp(cmap, ...)
}
