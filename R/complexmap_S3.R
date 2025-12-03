#' Internal Constructor for the ComplexMap S3 Object
#'
#' @description
#' Creates a `ComplexMap` object.
#'
#' @details
#' **Systems Biology Rationale:**
#' This constructor validates that the "Physical" (Nodes/Edges) and "Functional"
#' (Attributes) layers are correctly merged before finalizing the object.
#'
#' @param nodes A `tibble` containing node attributes and layout.
#' @param edges A `tibble` containing the network edge list.
#'
#' @return A `ComplexMap` S3 object.
#'
#' @keywords internal
#' @noRd
.new_ComplexMap <- function(nodes, edges) {
  # 1. Validation
  required_node_cols <- c("complexId", "primaryFunctionalDomain", "colorHex")
  missing_cols <- setdiff(required_node_cols, names(nodes))
  
  if (length(missing_cols) > 0) {
    stop("Invalid ComplexMap node table. Missing columns: ", 
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  # 2. Construction
  structure(
    list(
      nodes = nodes,
      edges = edges
    ),
    class = "ComplexMap"
  )
}

#' Print Method for a ComplexMap Object
#'
#' @description
#' Provides a "Systems Biology Dashboard" summary of the ComplexMap object.
#'
#' @details
#' This summary is designed to give immediate feedback on the **Functional Diversity**
#' of the landscape:
#' 
#' - **Functional Diversity:** The count of distinct functional labels (colors). 
#'   A high number (e.g., >15) indicates a rich, specific landscape. A low number
#'   indicates over-clustering or generic annotation.
#' - **Annotation Coverage:** The percentage of complexes that could be assigned
#'   a function. High coverage with high diversity is the ideal state.
#'
#' @param x A `ComplexMap` object.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the object.
#'
#' @export
print.ComplexMap <- function(x, ...) {
  # Extract Data
  numNodes <- nrow(x$nodes)
  numEdges <- nrow(x$edges)
  
  # Diversity Metrics
  valid_domains <- x$nodes$primaryFunctionalDomain[
    x$nodes$primaryFunctionalDomain != "Unenriched" & 
      !is.na(x$nodes$primaryFunctionalDomain)
  ]
  
  numThemes <- length(unique(valid_domains))
  coverage_pct <- if (numNodes > 0) (length(valid_domains) / numNodes) * 100 else 0
  
  # Physical Density
  density_str <- if (numNodes > 0) {
    sprintf("%.2f edges/node", numEdges / numNodes)
  } else "N/A"
  
  # Dashboard Output
  cat("# ComplexMap Object (Physical-First Layout)\n")
  cat(sprintf("# \u2500\u2500 Physical Structure: %d nodes, %d edges (%s)\n", 
              numNodes, numEdges, density_str))
  
  cat(sprintf("# \u2500\u2500 Functional Landscape:\n"))
  cat(sprintf("#    \u2022 Diversity: %d distinct functional domains (colors)\n", numThemes))
  cat(sprintf("#    \u2022 Coverage:  %.1f%% of complexes annotated\n", coverage_pct))
  
  # Hints
  cat("# \u2500\u2500 Accessors: `getNodeTable()`, `getEdgeTable()`\n")
  cat("# \u2500\u2500 Analysis:  `summarizeThemes()` to identify physical machines.\n")
  
  invisible(x)
}

#' Get the Node Table from a ComplexMap Object
#'
#' @description
#' Access the node attributes (layout, function, specificity scores).
#'
#' @param cm A `ComplexMap` object.
#' @return A `tibble`.
#' @export
getNodeTable <- function(cm) {
  return(cm$nodes)
}

#' Get the Edge Table from a ComplexMap Object
#'
#' @description
#' Access the network edge list (weights, similarity types).
#'
#' @param cm A `ComplexMap` object.
#' @return A `tibble`.
#' @export
getEdgeTable <- function(cm) {
  return(cm$edges)
}