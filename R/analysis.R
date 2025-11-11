utils::globalVariables(c("themeId", "primaryFunctionalDomain"))
#' Summarize Major Biological Themes in a Complex Map
#'
#' @description
#' This function analyzes the network structure within a `ComplexMap` object to
#' identify and summarize major biological themes. It uses community detection
#' algorithms to find densely connected clusters of nodes (themes).
#'
#' @details
#' The function performs the following steps:
#' 1.  It constructs an `igraph` graph object from the node and edge tables of
#'     the `complexMapObject`.
#' 2.  It applies a community detection algorithm (e.g., Louvain, as default)
#'     to partition the network into clusters or "themes".
#' 3.  For each identified theme, it generates a descriptive `themeLabel` by
#'     finding the most frequently occurring `primaryFunctionalDomain` among
#'     the member nodes (excluding "Unenriched").
#' 4.  It calculates summary statistics for each theme, including the number of
#'     nodes and edges it contains.
#'
#' @param complexMapObject A `ComplexMap` object returned by
#'   `createComplexMap()`.
#' @param method A character string specifying the community detection algorithm
#'   to use. Must be a valid `igraph` clustering function (e.g., "louvain",
#'   "walktrap", "infomap"). Defaults to "louvain".
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A `tibble` where each row represents a summarized theme. The tibble
#'   contains the following columns:
#'   - `themeId`: A unique integer identifier for the theme.
#'   - `themeLabel`: A descriptive label for the theme.
#'   - `nodeCount`: The number of nodes (complexes) in the theme.
#'   - `edgeCount`: The number of internal edges within the theme.
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
summarizeThemes <- function(complexMapObject, method="louvain", verbose=TRUE) {
  if (verbose) {
    message(
      sprintf("Summarizing themes using the '%s' community algorithm...",
              method)
    )
  }
  
  if (!requireNamespace("igraph", quietly=TRUE)) {
    stop("Package 'igraph' is required for theme summarization.", call.=FALSE)
  }
  
  nodes <- getNodeTable(complexMapObject)
  edges <- getEdgeTable(complexMapObject)
  
  if (nrow(nodes) == 0 || nrow(edges) == 0) {
    warning("Cannot summarize themes with empty node or edge lists.")
    return(tibble::tibble(
      themeId=integer(),
      themeLabel=character(),
      nodeCount=integer(),
      edgeCount=integer()
    ))
  }
  
  # 1. Create an igraph object
  # Ensure node IDs are characters to match edge list
  nodes$complexId <- as.character(nodes$complexId)
  graph <- igraph::graph_from_data_frame(
    d=edges, vertices=nodes, directed=FALSE
  )
  
  # 2. Perform community detection using the specified method
  communityAlgorithm <- switch(
    method,
    "louvain" = igraph::cluster_louvain,
    "walktrap" = igraph::cluster_walktrap,
    "infomap" = igraph::cluster_infomap,
    stop("Unsupported community detection method specified.")
  )
  communities <- communityAlgorithm(graph)
  nodes$themeId <- igraph::membership(communities)
  
  # 3. For each community, generate a theme label and summary stats
  themeSummaries <- nodes %>%
    dplyr::group_by(themeId) %>%
    dplyr::summarise(
      # Find the most common functional domain for the label
      themeLabel = {
        validDomains <- primaryFunctionalDomain[
          primaryFunctionalDomain != "Unenriched"
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
  
  # 4. Calculate the internal edge count for each theme
  subgraphEdgeCounts <- vapply(
    seq_along(communities),
    function(id) {
      subgraph <- igraph::induced_subgraph(
        graph, which(igraph::membership(communities) == id)
      )
      as.integer(igraph::ecount(subgraph))
    },
    integer(1)
  )
  themeSummaries$edgeCount <- subgraphEdgeCounts[themeSummaries$themeId]
  
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
#' - `type = "protein"`: The `query` is a character string representing a
#'   protein ID. The function searches the `proteins` column in the node table
#'   and returns all complexes that contain that protein. The search is exact
#'   and case-sensitive.
#' - `type = "complex"`: The `query` is a character string for a complex ID
#'   (e.g., "CpxMap_0001"). It returns the specific row from the node table
#'   corresponding to that complex.
#' - `type = "theme"`: The `query` is a character string corresponding to a
#'   `themeLabel` from `summarizeThemes()`. To use this, the themes must first
#'   be calculated and their IDs added to the node table. If the `themeId`
#'   column is not found, the function will stop with an informative error.
#'
#' @param complexMapObject A `ComplexMap` object returned by
#'   `createComplexMap()`.
#' @param query A character string containing the search term (a protein ID,
#'   complex ID, or theme label).
#' @param type A character string specifying the type of query. Must be one of
#'   `"protein"`, `"complex"`, or `"theme"`.
#'
#' @return A `tibble` containing the rows from the node table that match the
#'   query. Returns an empty tibble if no matches are found.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Assume 'cm_obj' is a valid ComplexMap object
#'
#' # Query for a specific protein (e.g., "UBA1")
#' # uba1_complexes <- queryMap(cm_obj, query = "UBA1", type = "protein")
#'
#' # Query for a specific complex ID
#' # cpx1_data <- queryMap(cm_obj, query = "CpxMap_0001", type = "complex")
#'
#' # To query by theme, you must first run summarizeThemes and add the IDs
#' # theme_summary <- summarizeThemes(cm_obj)
#' # nodes <- getNodeTable(cm_obj)
#' # graph <- igraph::graph_from_data_frame(getEdgeTable(cm_obj))
#' # communities <- igraph::cluster_louvain(graph)$membership
#' # nodes$themeId <- communities
#' # cm_obj$nodes <- nodes # Update the object
#' #
#' # Then you could query by a theme label from the summary
#' # proteasome_nodes <- queryMap(cm_obj, query = "Proteasome", type = "theme")
#'
queryMap <- function(complexMapObject, query, type) {
  nodes <- getNodeTable(complexMapObject)
  
  result <- switch(
    type,
    "protein" = {
      # Use stringr for robust pattern matching. We look for the query as a
      # whole word to avoid matching substrings (e.g., "UBA1" shouldn't
      # match "UBA10"). The pattern `\\b` matches a word boundary.
      pattern <- paste0("\\b", query, "\\b")
      nodes[stringr::str_detect(nodes$proteins, pattern), ]
    },
    
    "complex" = {
      nodes[nodes$complexId == query, ]
    },
    
    "theme" = {
      # This query type requires theme analysis to have been run first.
      if (!"themeLabel" %in% names(nodes) && !"themeId" %in% names(nodes)) {
        stop("Cannot query by theme. Please run `summarizeThemes()` first ",
             "and add theme information to the ComplexMap object.",
             call.=FALSE)
      }
      # This assumes a 'themeLabel' column exists or can be derived.
      # For simplicity, we require a column named 'themeLabel' for now.
      if (!"themeLabel" %in% names(nodes)) {
        stop("Node table is missing the 'themeLabel' column required for this query.",
             call.=FALSE)
      }
      nodes[nodes$themeLabel == query, ]
    },
    
    # Default case for invalid 'type'
    stop(sprintf("Invalid query type '%s'. Must be one of 'protein', ",
                 "'complex', or 'theme'.", type), call.=FALSE)
  )
  
  # Ensure a consistent return type (tibble) even if there are no matches
  if (nrow(result) == 0) {
    warning("Query returned no results for '", query, "'.")
    return(tibble::tibble())
  }
  
  return(result)
}