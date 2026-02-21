utils::globalVariables(c("primaryFunctionalDomain", "nodeCount", "purity"))

#' Query a ComplexMap Object
#'
#' @description
#' Targeted querying of complexes or proteins.
#'
#' @details
#' - `type="protein"`: Regex search for protein members.
#' - `type="complex"`: Exact match for Complex ID.
#'
#' @param complexMapObject A `ComplexMap` object.
#' @param query Search string.
#' @param type "protein" or "complex".
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
    stop("Invalid query type. Use 'protein' or 'complex'.")
  )
  
  if (nrow(result) == 0) warning(sprintf("No matches for '%s'.", query), call. = FALSE)
  return(result)
}