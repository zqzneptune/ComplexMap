#' Export a ComplexMap Network for External Tools
#'
#' @description
#' Exports the `ComplexMap` nodes and edges to compatible file formats, preserving
#' systems biology attributes like Specificity Scores and Functional Domains.
#'
#' @details
#' **Systems Biology Rationale:**
#' To verify the landscape in external tools (e.g., Cytoscape Desktop), this function
#' ensures that the specific calculated attributes are correctly formatted:
#'
#' - **`score`**: The Specificity Score.
#' - **`primaryFunctionalDomain`**: The Specific Label.
#' - **`colorHex`**: The calculated blended color.
#'
#' It performs a "Sanitization" step to convert any R-specific list columns
#' into semi-colon separated strings and fills `NA` values to ensure safe import.
#'
#' @param complexMapObject A `ComplexMap` object.
#' @param filePrefix Character string for output filenames (e.g., "results/my_map").
#' @param format Output format. One of:
#'   - `"cytoscape"`: Writes separate node and edge TSV files for Cytoscape Desktop.
#'   - `"graphml"`: Writes a GraphML file (via `igraph`).
#'   - `"tsv"`: Writes a single node-attribute TSV file.
#' @param verbose Logical.
#'
#' @return None (Writes files to disk).
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
exportNetwork <- function(complexMapObject, filePrefix, format = "cytoscape",
                          verbose = TRUE) {
  
  # 1. Extract Data
  nodes <- getNodeTable(complexMapObject)
  edges <- getEdgeTable(complexMapObject)
  
  if (nrow(nodes) == 0) {
    warning("ComplexMap object is empty. Nothing to export.")
    return(invisible(NULL))
  }
  
  # 2. Sanitize Data for Export
  # Helper to flatten list columns and handle NAs
  sanitize_df <- function(df) {
    df %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(is.list), 
        ~ vapply(., paste, collapse = "; ", FUN.VALUE = character(1))
      )) %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(is.character), 
        ~ tidyr::replace_na(., "")
      )) %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(is.numeric),
        ~ tidyr::replace_na(., 0)
      ))
  }
  
  nodes_clean <- sanitize_df(nodes)
  edges_clean <- sanitize_df(edges)
  
  if (format == "cytoscape") {
    nodeFile <- paste0(filePrefix, "_nodes.tsv")
    edgeFile <- paste0(filePrefix, "_edges.tsv")

    if (verbose) {
      message("Exporting to Cytoscape format (node + edge TSV)...")
      message(sprintf(" -> Nodes: %s (%d records)", nodeFile, nrow(nodes_clean)))
      message(sprintf(" -> Edges: %s (%d records)", edgeFile, nrow(edges_clean)))
    }

    utils::write.table(nodes_clean, file = nodeFile, sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(edges_clean, file = edgeFile, sep = "\t", row.names = FALSE, quote = FALSE)

  } else if (format == "graphml") {
    outFile <- paste0(filePrefix, ".graphml")

    if (verbose) message(sprintf("Exporting to GraphML: %s", outFile))

    graphObj <- igraph::graph_from_data_frame(
      d = edges_clean[
        c("source_complex_id", "target_complex_id",
          setdiff(names(edges_clean), c("source_complex_id", "target_complex_id")))],
      vertices = nodes_clean,
      directed = FALSE
    )
    igraph::write_graph(graphObj, file = outFile, format = "graphml")

  } else if (format == "tsv") {
    nodeFile <- paste0(filePrefix, "_nodes.tsv")

    if (verbose) message(sprintf("Exporting node TSV: %s (%d records)", nodeFile, nrow(nodes_clean)))
    utils::write.table(nodes_clean, file = nodeFile, sep = "\t", row.names = FALSE, quote = FALSE)

  } else {
    stop(sprintf("Unsupported format '%s'. Use 'cytoscape', 'graphml', or 'tsv'.", format))
  }

  if (verbose) message("Export complete.")
}