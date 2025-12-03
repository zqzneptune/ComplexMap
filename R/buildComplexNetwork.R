utils::globalVariables(c("compSim", "funcSim"))

#' Calculate a Similarity Score Between Two Sets
#'
#' @description
#' Internal helper to compute Jaccard, Overlap, or Dice coefficients.
#'
#' @keywords internal
.calculateSimilarity <- function(set1, set2, method="jaccard") {
  len1 <- length(set1)
  len2 <- length(set2)
  
  if (len1 == 0 || len2 == 0) return(0)
  
  intersectionSize <- length(intersect(set1, set2))
  if (intersectionSize == 0) return(0)
  
  # Note: "overlap" (Simpson) is available but discouraged for 
  # diverse landscape generation as it merges subsets into parents.
  score <- switch(
    method,
    "jaccard" = intersectionSize / (len1 + len2 - intersectionSize),
    "overlap" = intersectionSize / min(len1, len2),
    "dice"    = (2 * intersectionSize) / (len1 + len2),
    stop("Invalid similarity method. Use 'jaccard', 'overlap', or 'dice'.")
  )
  return(score)
}

#' Build a Complex-Complex Interaction Network
#'
#' @description
#' Constructs a network of complexes where edges represent similarity. 
#' 
#' @details
#' **Systems Biology Rationale:**
#' To generate a diverse and physically meaningful landscape, this function
#' defaults to prioritizing **Compositional Similarity** (shared proteins) over
#' functional similarity.
#' 
#' - **Composition (Protein Identity)** is treated as the "Physical Truth".
#' - **Function (Enrichment)** is treated as an attribute.
#' 
#' The `alpha` parameter controls this balance. The new default (0.75) gives
#' 75% weight to physical overlaps. This prevents generic functional terms
#' (like "Cell Cycle") from artificially pulling distinct physical complexes
#' into a single cluster, while still allowing functionally related complexes
#' to drift closer than unrelated ones.
#'
#' @param complexes A named list of protein complexes.
#' @param enrichments A named list of enrichment results (from `runComplexEnrichment`).
#' @param mode Edge weight mode: "compositional", "functional", or "combined".
#'   Defaults to "combined".
#' @param similarityMethod The metric for similarity. Defaults to **"jaccard"**.
#'   **Warning:** Using "overlap" is discouraged as it tends to collapse 
#'   diverse hierarchies into single blobs.
#' @param alpha Numeric (0-1). Weight given to Compositional Similarity in 
#'   "combined" mode. Defaults to **0.75** (favoring physical structure).
#' @param nCores Number of cores for parallel processing.
#' @param chunkSize Size of processing chunks.
#' @param verbose Logical.
#'
#' @return A `tibble` of network edges with `weight` representing the
#'   calculated similarity.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
buildComplexNetwork <- function(complexes, enrichments,
                                mode="combined",
                                similarityMethod="jaccard",
                                alpha=0.75,
                                nCores=NULL,
                                chunkSize=1000,
                                verbose=TRUE) {
  
  if (verbose) {
    message(sprintf("Building complex network (%s mode, alpha=%.2f)...", mode, alpha))
    if (similarityMethod == "overlap") {
      message("! NOTE: 'overlap' method selected. This may reduce landscape diversity.")
    }
  }
  
  if (length(complexes) < 2) {
    if (verbose) message("Fewer than 2 complexes; skipping network construction.")
    return(tibble::tibble(
      source_complex_id = character(), target_complex_id = character(),
      compSim = double(), funcSim = double(),
      sharedProt = integer(), sharedFunc = integer(),
      weight = double(), similarity_mode = character()
    ))
  }
  
  # --- 1. Setup Parallel Processing ---
  if (is.null(nCores)) nCores <- max(1, parallel::detectCores() - 1)
  
  # Future-proof parallel plan
  oldPlan <- future::plan(future::multisession, workers=nCores)
  on.exit(future::plan(oldPlan), add=TRUE)
  
  # --- 2. Prepare Data ---
  complexIds <- names(complexes)
  nComplexes <- length(complexIds)
  pairIndices <- utils::combn(nComplexes, 2, simplify=FALSE)
  
  # Split into chunks
  if (length(pairIndices) > 0) {
    nChunks <- ceiling(length(pairIndices) / chunkSize)
    groupFactor <- rep(seq_len(nChunks), each=chunkSize, length.out=length(pairIndices))
    pairChunks <- split(pairIndices, groupFactor)
  } else {
    pairChunks <- list()
  }
  
  if (verbose) message(sprintf("Processing %d pairs using %d cores...", length(pairIndices), nCores))
  
  # --- 3. Core Processing Function ---
  processChunk <- function(chunkPairs) {
    chunkResults <- vector("list", length(chunkPairs))
    
    for (idx in seq_along(chunkPairs)) {
      pairIdx <- chunkPairs[[idx]]
      i <- pairIdx[1]; j <- pairIdx[2]
      c1Id <- complexIds[i]; c2Id <- complexIds[j]
      
      # Compositional Similarity (Protein Overlap)
      # Always Jaccard-based logic implicitly via the chosen method
      compSim <- .calculateSimilarity(complexes[[c1Id]], complexes[[c2Id]], method=similarityMethod)
      
      # Functional Similarity (Enrichment Overlap)
      # Handle cases where complexes have NO enrichments gracefully
      terms1 <- if(is.null(enrichments[[c1Id]])) character(0) else enrichments[[c1Id]]$ID
      terms2 <- if(is.null(enrichments[[c2Id]])) character(0) else enrichments[[c2Id]]$ID
      
      funcSim <- .calculateSimilarity(terms1, terms2, method=similarityMethod)
      
      # Counts for metadata
      sharedProt <- length(intersect(complexes[[c1Id]], complexes[[c2Id]]))
      sharedFunc <- length(intersect(terms1, terms2))
      
      chunkResults[[idx]] <- list(
        source_complex_id = c1Id,
        target_complex_id = c2Id,
        compSim = compSim,
        funcSim = funcSim,
        sharedProt = sharedProt,
        sharedFunc = sharedFunc
      )
    }
    return(chunkResults)
  }
  
  # --- 4. Execution ---
  if (verbose && requireNamespace("progressr", quietly=TRUE)) {
    progressr::with_progress({
      p <- progressr::progressor(steps=length(pairChunks))
      chunkList <- future.apply::future_lapply(pairChunks, function(chunk) {
        res <- processChunk(chunk)
        p()
        res
      }, future.seed=TRUE)
    })
  } else {
    chunkList <- future.apply::future_lapply(pairChunks, processChunk, future.seed=TRUE)
  }
  
  # --- 5. Assembly ---
  if (length(chunkList) > 0) {
    allResults <- unlist(chunkList, recursive=FALSE)
    finalNetworkDf <- dplyr::bind_rows(allResults)
  } else {
    finalNetworkDf <- tibble::tibble(source_complex_id=character(), target_complex_id=character())
  }
  
  # --- 6. Weight Calculation ---
  finalNetworkDf <- finalNetworkDf %>%
    dplyr::mutate(
      weight = dplyr::case_when(
        mode == "compositional" ~ compSim,
        mode == "functional"   ~ funcSim,
        mode == "combined"     ~ (alpha * compSim) + ((1 - alpha) * funcSim),
        TRUE                   ~ NA_real_
      ),
      similarity_mode = mode
    ) %>%
    dplyr::filter(weight > 0)
  
  if (verbose) {
    message(sprintf("Network built: %d edges retained.", nrow(finalNetworkDf)))
  }
  
  return(finalNetworkDf)
}