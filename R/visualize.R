# visualize.R — DEPRECATED
# All static visualization functions have been removed.
# Use explore(cmap) to launch the interactive Cytoscape.js explorer.

#' @keywords internal
.visualize_deprecated <- function() {
  stop(
    "Static visualization functions have been removed in ComplexMap >= 2.0.\n",
    "Use `explore(cmap)` to launch the interactive Cytoscape.js explorer.",
    call. = FALSE
  )
}