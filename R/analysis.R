utils::globalVariables(c("themeId", "primaryFunctionalDomain", "themeLabel",
                         "nodeCount"))

#' Summarize Major Biological Themes in a Complex Map
#'
#' @description
#' This function analyzes the network structure within a `ComplexMap` object to
#' identify and summarize major biological themes. It uses community detection
#' algorithms to find densely connected clusters of nodes (themes).
#'
#' @details
#' The function performs the following steps:
#' 1.  It constructs an `igraph` graph object from the node and edge tables.
#' 
#' 2.  It applies a community detection algorithm (e.g., Louvain, as default)
#'     to partition the network into clusters or "themes".
#'     
#' 3.  For each theme, it generates a descriptive `themeLabel` by finding the
#'     most frequently occurring `primaryFunctionalDomain` among the member
#'     nodes (excluding "Unenriched").
#'     
#' 4.  It calculates summary statistics for each theme, including node and edge counts.
#'
#' @param complexMapObject A `ComplexMap` object returned by `createComplexMap()`.
#' @param method A character string specifying the community detection algorithm
#'   to use. Must be a valid `igraph` clustering function (e.g., "louvain",
#'   "walktrap", "infomap"). Defaults to "louvain".
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A `tibble` where each row represents a summarized theme, containing
#'   `themeId`, `themeLabel`, `nodeCount`, and `edgeCount`.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Assume 'cm_obj' is a valid ComplexMap object created by createComplexMap()
#' # if (requireNamespace("igraph", quietly = TRUE)) {
#' #   themeSummary <- summarizeThemes(cm_obj)
#' #   print(themeSummary)
#' # }
#'
summarizeThemes <- function(complexMapObject, method = "louvain", verbose = TRUE) {
  if (verbose) {
    message(
      sprintf(
        "Summarizing themes using the '%s' community algorithm...",
        method)
    )
  }
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for theme summarization.", call. = FALSE)
  }
  
  nodes <- getNodeTable(complexMapObject)
  edges <- getEdgeTable(complexMapObject)
  
  if (nrow(nodes) == 0 || nrow(edges) == 0) {
    warning("Cannot summarize themes with empty node or edge lists.")
    return(tibble::tibble(
      themeId = integer(),
      themeLabel = character(),
      nodeCount = integer(),
      edgeCount = integer()
    ))
  }
  
  nodes$complexId <- as.character(nodes$complexId)
  graph <- igraph::graph_from_data_frame(
    d = edges, vertices = nodes, directed = FALSE
  )
  
  communityAlgorithm <- switch(
    method,
    "louvain" = igraph::cluster_louvain,
    "walktrap" = igraph::cluster_walktrap,
    "infomap" = igraph::cluster_infomap,
    stop("Unsupported community detection method specified.")
  )
  communities <- communityAlgorithm(graph)
  
  # Safely create a mapping from complexId to themeId
  theme_mapping <- tibble::tibble(
    complexId = igraph::V(graph)$name,
    themeId = as.integer(igraph::membership(communities))
  )
  
  # Join mapping back to the original node table to get functional domains
  nodes_with_themes <- nodes %>%
    dplyr::left_join(theme_mapping, by = "complexId")
  
  themeSummaries <- nodes_with_themes %>%
    dplyr::filter(!is.na(themeId)) %>% # Exclude nodes not in the graph
    dplyr::group_by(themeId) %>%
    dplyr::summarise(
      themeLabel = {
        validDomains <- primaryFunctionalDomain[
          primaryFunctionalDomain != "Unenriched" & !is.na(primaryFunctionalDomain)
        ]
        if (length(validDomains) > 0) {
          names(which.max(table(validDomains)))[1]
        } else {
          "Mixed or Unenriched"
        }
      },
      nodeCount = dplyr::n(),
      .groups = "drop"
    )
  
  subgraphEdgeCounts <- vapply(
    unique(theme_mapping$themeId),
    function(id) {
      members <- theme_mapping$complexId[theme_mapping$themeId == id]
      subgraph <- igraph::induced_subgraph(graph, members)
      as.integer(igraph::ecount(subgraph))
    },
    integer(1)
  )
  
  edge_counts_df <- tibble::tibble(
    themeId = unique(theme_mapping$themeId),
    edgeCount = subgraphEdgeCounts
  )
  
  themeSummaries <- themeSummaries %>%
    dplyr::left_join(edge_counts_df, by = "themeId")
  
  if (verbose) {
    message(sprintf("Identified %d distinct themes.", nrow(themeSummaries)))
  }
  
  return(themeSummaries)
}


#' Query a ComplexMap Object for Specific Information
#'
#' @description
#' Allows for targeted querying of a `ComplexMap` object to find nodes
#' (complexes) that match specific criteria, such as containing a particular
#' protein or belonging to a biological theme.
#'
#' @details
#' This function supports three distinct modes of querying:
#' - `type = "protein"`: Searches for complexes containing a specific protein.
#' - `type = "complex"`: Retrieves a specific complex by its ID.
#' - `type = "theme"`: Finds all complexes belonging to a given theme. This
#'   requires theme information to be added to the `ComplexMap` object first
#'   (see the "Advanced Analysis" vignette for an example).
#'
#' @param complexMapObject A `ComplexMap` object.
#' @param query A character string containing the search term.
#' @param type A character string: one of `"protein"`, `"complex"`, or `"theme"`.
#'
#' @return A `tibble` containing the rows from the node table that match the
#'   query. Returns an empty tibble if no matches are found.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Assume 'cm_obj' is a valid ComplexMap object
#' # uba1_complexes <- queryMap(cm_obj, query = "UBA1", type = "protein")
#'
queryMap <- function(complexMapObject, query, type) {
  nodes <- getNodeTable(complexMapObject)
  
  result <- switch(
    type,
    "protein" = {
      pattern <- paste0("\\b", query, "\\b")
      nodes[stringr::str_detect(nodes$proteins, pattern), ]
    },
    
    "complex" = {
      nodes[nodes$complexId == query, ]
    },
    
    "theme" = {
      if (!"themeLabel" %in% names(nodes)) {
        stop("Cannot query by theme. Column 'themeLabel' not found in the node table.",
             call. = FALSE)
      }
      # --- ROBUST QUERY ---
      # Filter for the theme and also remove any rows where themeLabel is NA
      nodes[!is.na(nodes$themeLabel) & nodes$themeLabel == query, ]
    },
    
    stop(sprintf("Invalid query type '%s'. Must be one of 'protein', 'complex', or 'theme'.", type),
         call. = FALSE)
  )
  
  if (nrow(result) == 0) {
    warning("Query returned no results for '", query, "'.", call. = FALSE)
    return(tibble::tibble())
  }
  
  return(result)
}