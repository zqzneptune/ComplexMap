utils::globalVariables(c("themeId", "primaryFunctionalDomain", "themeLabel",
                         "nodeCount", "purity", "themePurity"))

#' Summarize Major Biological Themes in a Complex Map
#'
#' @description
#' Identifies and summarizes physical neighborhoods (themes) in the complex map
#' using community detection.
#'
#' @details
#' **Systems Biology Rationale:**
#' Since the network layout is driven primarily by physical composition (alpha=0.75),
#' the communities detected here represent **Physical Machines** or neighborhoods.
#' 
#' This function characterizes these machines by their functions. A single physical
#' machine might contain complexes with slightly different specific labels.
#' To reflect this, the labeling logic is:
#' 1.  Identify the most frequent functional domain in the cluster.
#' 2.  Calculate **Theme Purity** (what % of nodes share this domain).
#' 3.  If purity is low (< 50%), the label combines the top 2 domains
#'     (e.g., "Function A / Function B") to indicate multi-functionality.
#'
#' By default, this function adds the theme assignments directly to the node
#' table of the `ComplexMap` object, making it easy to visualize or query by theme.
#'
#' @param complexMapObject A `ComplexMap` object.
#' @param method Community detection method ("louvain", "walktrap", "infomap").
#'   Defaults to "louvain".
#' @param add_to_object Logical. If `TRUE` (default), returns the `ComplexMap`
#'   object with `themeId` and `themeLabel` columns added to the node table. If
#'   `FALSE`, returns a summary `tibble` of the themes.
#' @param verbose Logical.
#'
#' @return
#' If `add_to_object = TRUE`, a modified `ComplexMap` object.
#' If `add_to_object = FALSE`, a `tibble` with columns:
#'   - `themeId`: Cluster ID.
#'   - `themeLabel`: The consensus functional label.
#'   - `themePurity`: The fraction of nodes matching the primary label (0-1).
#'   - `nodeCount`: Number of complexes.
#'   - `edgeCount`: Number of internal edges.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
summarizeThemes <- function(complexMapObject, method = "louvain",
                            add_to_object = TRUE, verbose = TRUE) {
  if (verbose) message(sprintf("Summarizing physical themes using '%s'...", method))
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }
  
  nodes <- getNodeTable(complexMapObject)
  edges <- getEdgeTable(complexMapObject)
  
  if (nrow(nodes) == 0) {
    warning("Cannot summarize themes for an empty map.")
    return(if (add_to_object) complexMapObject else tibble::tibble())
  }
  
  # Build graph
  nodes$complexId <- as.character(nodes$complexId)
  graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # Community Detection
  communityAlg <- switch(
    method,
    "louvain" = igraph::cluster_louvain,
    "walktrap" = igraph::cluster_walktrap,
    "infomap" = igraph::cluster_infomap,
    stop("Invalid method.")
  )
  communities <- communityAlg(graph)
  
  # Map themes
  theme_mapping <- tibble::tibble(
    complexId = igraph::V(graph)$name,
    themeId = as.integer(igraph::membership(communities))
  )
  
  nodes_with_themes <- nodes %>% dplyr::left_join(theme_mapping, by = "complexId")
  
  # Summarize with Diversity Logic
  themeSummaries <- nodes_with_themes %>%
    dplyr::filter(!is.na(themeId)) %>%
    dplyr::group_by(themeId) %>%
    dplyr::summarise(
      {
        valid <- primaryFunctionalDomain[primaryFunctionalDomain != "Unenriched" & !is.na(primaryFunctionalDomain)]
        
        if (length(valid) == 0) {
          lbl <- "Unenriched"
          pur <- 0
        } else {
          counts <- sort(table(valid), decreasing = TRUE)
          lbl <- names(counts)[1]
          pur <- as.numeric(counts[1]) / length(valid)
          
          # If high diversity (low purity), append second label
          if (pur < 0.50 && length(counts) > 1) {
            lbl <- paste(lbl, "/", names(counts)[2])
          }
        }
        tibble::tibble(themeLabel = lbl, themePurity = round(pur, 2))
      },
      nodeCount = dplyr::n(),
      .groups = "drop"
    )
  
  # Add Edge Counts
  edge_counts <- vapply(themeSummaries$themeId, function(id) {
    mems <- theme_mapping$complexId[theme_mapping$themeId == id]
    if(length(mems) < 2) return(0L)
    as.integer(igraph::ecount(igraph::induced_subgraph(graph, mems)))
  }, integer(1))
  
  themeSummaries$edgeCount <- edge_counts
  
  if (verbose) message(sprintf("Identified %d themes. Avg Purity: %.2f", 
                               nrow(themeSummaries), mean(themeSummaries$themePurity, na.rm = TRUE)))
  
  if (add_to_object) {
    # Join theme labels back to the full theme mapping
    full_theme_info <- theme_mapping %>%
      dplyr::left_join(
        dplyr::select(themeSummaries, themeId, themeLabel, themePurity),
        by = "themeId"
      )
    
    # Update the node table in the object
    updated_nodes <- nodes %>%
      dplyr::left_join(full_theme_info, by = "complexId")
    
    complexMapObject$nodes <- updated_nodes
    return(complexMapObject)
  } else {
    return(themeSummaries)
  }
}


#' Query a ComplexMap Object
#'
#' @description
#' Targeted querying of complexes, proteins, or themes.
#'
#' @details
#' - `type="protein"`: Regex search for protein members.
#' - `type="complex"`: Exact match for Complex ID.
#' - `type="theme"`: Exact match for Theme Label. Requires `summarizeThemes()`
#'   to be run on the object first with `add_to_object = TRUE`.
#'
#' @param complexMapObject A `ComplexMap` object.
#' @param query Search string.
#' @param type "protein", "complex", or "theme".
#'
#' @return A tibble of matching nodes.
#'
#' @export
queryMap <- function(complexMapObject, query, type) {
  nodes <- getNodeTable(complexMapObject)
  
  result <- switch(
    type,
    "protein" = nodes[stringr::str_detect(nodes$proteins, paste0("\\b", query, "\\b")), ],
    "complex" = nodes[nodes$complexId == query, ],
    "theme"   = {
      if (!"themeLabel" %in% names(nodes)) {
        stop("Cannot query by 'theme'. Please run `summarizeThemes()` on the object first.", call. = FALSE)
      }
      # Robust check for NA
      nodes[!is.na(nodes$themeLabel) & nodes$themeLabel == query, ]
    },
    stop("Invalid query type. Use 'protein', 'complex', or 'theme'.")
  )
  
  if (nrow(result) == 0) warning(sprintf("No matches for '%s'.", query), call. = FALSE)
  return(result)
}