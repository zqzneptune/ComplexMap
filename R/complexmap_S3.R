#' Internal Constructor for the ComplexMap S3 Object
#'
#' @description
#' Creates a new `ComplexMap` object from nodes and edges tibbles. This is a
#' low-level constructor and does not perform any validation. It is intended
#' for internal use only.
#'
#' @param nodes A `tibble` containing final node attributes and layout data from
#'   `computeMapTopology`.
#' @param edges A `tibble` containing the network edge list from
#'   `buildComplexNetwork`.
#'
#' @return A `ComplexMap` S3 object, which is a list containing the `nodes`
#'   and `edges` tibbles.
#'
#' @keywords internal
#' @noRd
.new_ComplexMap <- function(nodes, edges) {
  # Per the brief, the S3 object is a list with two named elements.
  # The structure() function is a clean and standard way to assign a class
  # attribute while creating the object.
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
#' Provides a user-friendly summary of the `ComplexMap` object when it is
#' printed to the console. This method is automatically called by the `print`
#' generic function.
#'
#' @param x A `ComplexMap` object.
#' @param ... Additional arguments passed to `print` (not used by this method).
#'
#' @return Invisibly returns the original `ComplexMap` object, allowing for
#'   its use in pipelines.
#'
#' @export
print.ComplexMap <- function(x, ...) {
  # Extract the number of nodes and edges directly from the object's tibbles.
  numNodes <- nrow(x$nodes)
  numEdges <- nrow(x$edges)
  
  # To count themes, we find the number of unique functional domains,
  # making sure to exclude the generic "Unenriched" category.
  # The `primaryFunctionalDomain` column is required for this.
  if ("primaryFunctionalDomain" %in% names(x$nodes)) {
    numThemes <- length(unique(
      x$nodes$primaryFunctionalDomain[
        x$nodes$primaryFunctionalDomain != "Unenriched"
      ]
    ))
  } else {
    # Provide a fallback if the column is somehow missing.
    numThemes <- 0
  }
  
  # Use cat() to print the formatted summary string to the console.
  # The "──" character (Unicode U+2500 U+2500) is used for styling.
  cat("# A ComplexMap Object\n")
  cat(sprintf("# \u2500\u2500 %d nodes and %d edges\n", numNodes, numEdges))
  cat(sprintf(
    "# \u2500\u2500 %d major biological themes identified.\n", numThemes
  ))
  cat("# \u2500\u2500 Use `getNodeTable()` or `getEdgeTable()` to access data.\n")
  
  # It is standard practice for print methods to return the original object
  # invisibly, which prevents it from being printed twice in a row.
  invisible(x)
}

#' Get the Node Table from a ComplexMap Object
#'
#' @description
#' A simple accessor function to extract the tibble of node attributes and
#' layout data from a `ComplexMap` object.
#'
#' @param cm A `ComplexMap` object, typically the output of
#'   `createComplexMap()`.
#'
#' @return A `tibble` containing the node data.
#'
#' @export
#' @examples
#' # Assuming 'myComplexMap' is a valid ComplexMap object
#' # nodeData <- getNodeTable(myComplexMap)
getNodeTable <- function(cm) {
  # This directly accesses the 'nodes' element of the list.
  # No validation is needed here, as it assumes a valid ComplexMap object.
  return(cm$nodes)
}

#' Get the Edge Table from a ComplexMap Object
#'
#' @description
#' A simple accessor function to extract the tibble of network edges from a
#' `ComplexMap` object.
#'
#' @param cm A `ComplexMap` object, typically the output of
#'   `createComplexMap()`.
#'
#' @return A `tibble` containing the edge data.
#'
#' @export
#' @examples
#' # Assuming 'myComplexMap' is a valid ComplexMap object
#' # edgeData <- getEdgeTable(myComplexMap)
getEdgeTable <- function(cm) {
  # This directly accesses the 'edges' element of the list.
  return(cm$edges)
}