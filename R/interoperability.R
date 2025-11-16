#' Export a ComplexMap Network for External Tools
#'
#' @description
#' This function exports the node and edge tables from a `ComplexMap` object
#' into a format suitable for widely used network visualization software, such
#' as Cytoscape.
#'
#' @details
#' The function currently supports one format:
#' 
#' - `"cytoscape"`: This option generates two separate tab-separated value
#'   (.tsv) files. 
#'   
#'   One file contains the node attributes (`<filePrefix>_nodes.tsv`)
#'   and the other contains the edge list with its attributes
#'   (`<filePrefix>_edges.tsv`). These files can be directly imported into
#'   Cytoscape's network and attribute tables.
#'
#' The function uses `utils::write.table` for robust and standard-compliant
#' file writing.
#'
#' @param complexMapObject A `ComplexMap` object returned by
#'   `createComplexMap()`.
#' @param filePrefix A character string that will be used as the prefix for the
#'   output filenames. For example, a `filePrefix` of "my_map" will result in
#'   "my_map_nodes.tsv" and "my_map_edges.tsv".
#' @param format A character string specifying the output format. Currently,
#'   only `"cytoscape"` is supported.
#' @param verbose A logical value indicating whether to print a confirmation
#'   message upon successful export.
#'
#' @return The function is called for its side effect of writing files to disk.
#'   It does not return a value.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Assume 'cm_obj' is a valid ComplexMap object
#' # dir <- tempdir() # Use a temporary directory for the example
#' # exportNetwork(cm_obj, filePrefix = file.path(dir, "myNetwork"))
#'
exportNetwork <- function(complexMapObject, filePrefix, format="cytoscape",
                          verbose=TRUE) {
  # 1. Extract the node and edge tibbles from the object
  nodes <- getNodeTable(complexMapObject)
  edges <- getEdgeTable(complexMapObject)
  
  if (format == "cytoscape") {
    nodeFile <- paste0(filePrefix, "_nodes.tsv")
    edgeFile <- paste0(filePrefix, "_edges.tsv")
    
    if (verbose) {
      message(
        "Exporting network in Cytoscape format to:\n -> ",
        nodeFile, "\n -> ", edgeFile
      )
    }
    
    # 2. Write the nodes table using utils::write.table
    utils::write.table(
      nodes,
      file = nodeFile,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    # 3. Write the edges table using utils::write.table
    utils::write.table(
      edges,
      file = edgeFile,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
  } else {
    stop(sprintf("Unsupported format '%s' specified.", format))
  }
  
  if (verbose) message("Export complete.")
}