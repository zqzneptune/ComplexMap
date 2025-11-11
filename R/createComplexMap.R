#' Create a Complete Complex Map Object
#'
#' @description
#' A high-level wrapper function that executes the entire `ComplexMap` workflow.
#' It takes a list of protein complexes and a gene set matrix (GMT) and performs
#' refinement, enrichment, network construction, and topology calculation.
#'
#' @details
#' This function serves as the primary entry point for most analyses. It
#' internally calls the core workflow in the following order:
#' 1.  `refineComplexList()`
#' 2.  `runComplexEnrichment()`
#' 3.  `buildComplexNetwork()`
#' 4.  `generateNodeAttributes()`
#' 5.  `computeMapTopology()`
#'
#' Arguments for the underlying functions can be passed directly to this
#' wrapper via the `...` parameter.
#'
#' @param complexList A named list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param gmt A named list where each element is a character vector of genes,
#'   representing a functional gene set (e.g., from a GMT file).
#' @param verbose A logical value indicating whether to print progress messages
#'   for the entire workflow. Defaults to `TRUE`.
#' @param ... Additional arguments to be passed down to the core functions.
#'   Common arguments include:
#'   - `minSize`, `maxSize`, `mergeThreshold` (for `refineComplexList`)
#'   - `pAdjustMethod`, `pValueCutoff` (for `runComplexEnrichment`)
#'   - `mode`, `similarityMethod`, `alpha` (for `buildComplexNetwork`)
#'
#' @return A validated `ComplexMap` S3 object containing the final node and
#'   edge tables.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Assume 'demoComplexes' and a 'gmt' object are loaded
#' # gmtPath <- getExampleGmt()
#' # gmt <- getGmtFromFile(gmtPath, verbose = FALSE)
#'
#' # Run the full workflow with custom parameters
#' # complexMapObject <- createComplexMap(
#' #   demoComplexes,
#' #   gmt,
#' #   verbose = TRUE,
#' #   minSize = 5,
#' #   mergeThreshold = 0.8,
#' #   pValueCutoff = 0.01,
#' #   mode = "combined"
#' # )
#' # print(complexMapObject)
#'
createComplexMap <- function(complexList, gmt, verbose=TRUE, ...) {
  # Capture all additional arguments passed via the ellipsis (...)
  dotArgs <- list(...)
  
  if (verbose) message("--- Starting ComplexMap Workflow ---")
  
  # --- 1. Refine Complex List ---
  if (verbose) message("\nStep 1: Refining complex list...")
  refineArgs <- list(complexList=complexList, verbose=verbose)
  # Append relevant arguments from ... if they exist
  refineParams <- c("minSize", "maxSize", "mergeThreshold")
  refineArgs <- c(refineArgs, dotArgs[names(dotArgs) %in% refineParams])
  refinedList <- do.call(refineComplexList, refineArgs)
  
  # --- 2. Run Enrichment Analysis ---
  if (verbose) message("\nStep 2: Running enrichment analysis...")
  enrichArgs <- list(complexList=refinedList, gmt=gmt, verbose=verbose)
  enrichParams <- c("pAdjustMethod", "pValueCutoff")
  enrichArgs <- c(enrichArgs, dotArgs[names(dotArgs) %in% enrichParams])
  enrichmentResults <- do.call(runComplexEnrichment, enrichArgs)
  
  # --- 3. Build Complex Network ---
  if (verbose) message("\nStep 3: Building complex network...")
  networkArgs <- list(
    complexes=refinedList,
    enrichments=enrichmentResults,
    verbose=verbose
  )
  networkParams <- c("mode", "similarityMethod", "alpha", "nCores", "chunkSize")
  networkArgs <- c(networkArgs, dotArgs[names(dotArgs) %in% networkParams])
  networkEdges <- do.call(buildComplexNetwork, networkArgs)
  
  # --- 4. Generate Node Attributes ---
  if (verbose) message("\nStep 4: Generating node attributes...")
  nodeAttrArgs <- list(
    complexes=refinedList,
    enrichments=enrichmentResults,
    verbose=verbose
  )
  # Note: similarityMethod can be used by both buildComplexNetwork and this
  nodeAttrParams <- "similarityMethod"
  nodeAttrArgs <- c(nodeAttrArgs, dotArgs[names(dotArgs) %in% nodeAttrParams])
  nodeAttributes <- do.call(generateNodeAttributes, nodeAttrArgs)
  
  # --- 5. Compute Map Topology ---
  if (verbose) message("\nStep 5: Computing map topology...")
  topoArgs <- list(
    nodeAttributes=nodeAttributes,
    network=networkEdges,
    verbose=verbose
  )
  finalNodes <- do.call(computeMapTopology, topoArgs)
  
  if (verbose) message("\n--- ComplexMap Workflow Complete ---")
  
  # --- 6. Construct and Return Final S3 Object ---
  .new_ComplexMap(nodes=finalNodes, edges=networkEdges)
}