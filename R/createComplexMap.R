#' Create a Complete Complex Map Object
#'
#' @description
#' A high-level wrapper function that executes the entire `ComplexMap` workflow.
#' It takes a list of protein complexes and a gene set matrix (GMT) and performs
#' refinement, enrichment, network construction, and topology calculation.
#'
#' @details
#' This function is tuned to generate a **functionally diverse landscape**.
#' It enforces the following logic:
#' 1.  **Refinement:** Uses Jaccard similarity to merge only highly redundant complexes,
#'     preserving biological variants (subsets/supersets).
#' 2.  **Enrichment:** Calculates 'Fold Enrichment' to prioritize specific biological
#'     functions over generic ones.
#' 3.  **Network:** Builds the layout primarily based on **Physical Composition**
#'     (alpha = 0.75), using functional annotations only to group related clusters,
#'     not to collapse them.
#'
#' @param complexList A named list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param gmt A named list where each element is a character vector of genes.
#' @param similarityMethod The metric used for comparing complexes and functional terms.
#'   Defaults to **"jaccard"** to penalize size differences and maintain diversity.
#'   Avoid "overlap" if you wish to see specific pathways distinct from generic ones.
#' @param alpha Numeric (0-1). The weight given to physical protein composition
#'   versus functional similarity in the network layout. Defaults to **0.75**
#'   (75\% physical, 25\% functional).
#' @param verbose A logical value indicating whether to print progress messages.
#' @param ... Additional arguments passed to core functions:
#'   - `minSize`, `maxSize`, `mergeThreshold` (for `refineComplexList`)
#'   - `pAdjustMethod`, `pValueCutoff` (for `runComplexEnrichment`)
#'   - `geneSetDb` (for `generateNodeAttributes` semantic clustering)
#'
#' @return A validated `ComplexMap` S3 object.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
createComplexMap <- function(complexList, gmt, 
                             similarityMethod = "jaccard",
                             alpha = 0.75,
                             verbose = TRUE, ...) {
  
  # Capture additional arguments
  dotArgs <- list(...)
  
  if (verbose) {
    message("--- Starting ComplexMap Workflow ---")
    message(sprintf("Parameters: similarity='%s', alpha=%.2f (Diversity Priority)", 
                    similarityMethod, alpha))
  }
  
  # --- 1. Refine Complex List ---
  if (verbose) message("\nStep 1: Refining complex list...")
  
  # Explicitly pass similarityMethod to ensure we don't accidentally use 
  # "matching_score" or "overlap" unless the user specifically requested it via overrides.
  refineArgs <- list(
    complexList = complexList, 
    similarityMethod = similarityMethod, # Enforce Jaccard by default
    verbose = verbose
  )
  
  # Allow overrides from ... (e.g., mergeThreshold)
  refineParams <- c("minSize", "maxSize", "mergeThreshold")
  refineArgs <- c(refineArgs, dotArgs[names(dotArgs) %in% refineParams])
  
  refinement_output <- do.call(refineComplexList, refineArgs)
  refinedList <- refinement_output$refinedComplexes
  
  # --- 2. Run Enrichment Analysis ---
  if (verbose) message("\nStep 2: Running enrichment analysis...")
  enrichArgs <- list(complexList = refinedList, gmt = gmt, verbose = verbose)
  enrichParams <- c("pAdjustMethod", "pValueCutoff")
  enrichArgs <- c(enrichArgs, dotArgs[names(dotArgs) %in% enrichParams])
  enrichmentResults <- do.call(runComplexEnrichment, enrichArgs)
  
  # --- 3. Build Complex Network ---
  if (verbose) message("\nStep 3: Building complex network...")
  networkArgs <- list(
    complexes = refinedList,
    enrichments = enrichmentResults,
    similarityMethod = similarityMethod, # Enforce Jaccard consistency
    alpha = alpha,                       # Enforce Physical priority
    verbose = verbose
  )
  networkParams <- c("mode", "nCores", "chunkSize")
  networkArgs <- c(networkArgs, dotArgs[names(dotArgs) %in% networkParams])
  networkEdges <- do.call(buildComplexNetwork, networkArgs)
  
  # --- 4. Generate Node Attributes ---
  if (verbose) message("\nStep 4: Generating node attributes...")
  nodeAttrArgs <- list(
    complexes = refinedList,
    enrichments = enrichmentResults,
    similarityMethod = similarityMethod, # Enforce Jaccard for color clustering
    verbose = verbose
  )
  # Pass geneSetDb if it exists in ... for semantic clustering
  nodeAttrParams <- c("geneSetDb")
  nodeAttrArgs <- c(nodeAttrArgs, dotArgs[names(dotArgs) %in% nodeAttrParams])
  nodeAttributes <- do.call(generateNodeAttributes, nodeAttrArgs)
  
  # --- 5. Compute Map Topology ---
  if (verbose) message("\nStep 5: Computing map topology...")
  topoArgs <- list(
    nodeAttributes = nodeAttributes,
    network = networkEdges,
    verbose = verbose
  )
  finalNodes <- do.call(computeMapTopology, topoArgs)
  
  if (verbose) message("\n--- ComplexMap Workflow Complete ---")
  
  # --- 6. Construct and Return Final S3 Object ---
  .new_ComplexMap(nodes = finalNodes, edges = networkEdges)
}