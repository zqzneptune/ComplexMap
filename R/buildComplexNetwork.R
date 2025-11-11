utils::globalVariables(c("compSim", "funcSim"))
#' Calculate a Similarity Score Between Two Sets
#'
#' @description
#' An internal helper function that computes one of three common similarity
#' metrics: Jaccard, Overlap, or Dice.
#'
#' @param set1 A character vector representing the first set.
#' @param set2 A character vector representing the second set.
#' @param method The similarity metric to use. One of "jaccard", "overlap",
#'   or "dice".
#'
#' @return A numeric similarity score between 0 and 1. Returns 0 if either
#'   set is empty or if there is no intersection.
#'
#' @keywords internal
#'
.calculateSimilarity <- function(set1, set2, method="jaccard") {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  
  intersectionSize <- length(intersect(set1, set2))
  if (intersectionSize == 0) {
    return(0)
  }
  
  score <- switch(
    method,
    "jaccard" = {
      intersectionSize / length(union(set1, set2))
    },
    "overlap" = {
      intersectionSize / min(length(set1), length(set2))
    },
    "dice" = {
      (2 * intersectionSize) / (length(set1) + length(set2))
    },
    stop("Invalid similarity method specified.")
  )
  return(score)
}

#' Build a Complex-Complex Interaction Network
#'
#' @description
#' Constructs a network of complexes where edges represent similarity. The
#' similarity can be based on shared proteins (compositional), shared
#' functional annotations (functional), or a weighted combination of both.
#'
#' @details
#' This function calculates all pairwise similarities between complexes in a
#' list. It is highly optimized for large datasets by using chunking and
#' parallel processing via the `future` framework.
#'
#' The final edge weight is determined by the `mode`:
#' - `"compositional"`: Uses only the protein similarity score.
#' - `"functional"`: Uses only the functional annotation similarity score.
#' - `"combined"`: Uses a weighted average:
#'   `alpha * compositional + (1 - alpha) * functional`.
#'
#' If the `progressr` package is installed, a progress bar will be displayed
#' during the parallel computation when `verbose = TRUE`.
#'
#' @param complexes A named list of protein complexes.
#' @param enrichments A named list of enrichment results, corresponding to the
#'   `complexes`. Typically the output of `runComplexEnrichment`.
#' @param mode A character string specifying how to calculate the final edge
#'   weight. One of "functional", "compositional", or "combined".
#' @param similarityMethod The metric used for both compositional and functional
#'   similarity. One of "jaccard", "overlap", or "dice".
#' @param alpha A numeric value (0-1) used only in "combined" mode to weigh
#'   the compositional similarity score.
#' @param nCores The number of CPU cores for parallel processing. Defaults to
#'   one less than available.
#' @param chunkSize The number of complex pairs to process in each parallel
#'   chunk.
#' @param verbose A logical value indicating whether to show progress messages
#'   and a progress bar.
#'
#' @return
#' A `tibble` representing the network edges. Each row includes the source
#' and target_complex_id complexes, their similarity scores, shared component
#' counts, the final calculated `weight`, and the `similarity_mode`.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data (from previous examples) ---
#' complexes <- list(
#'   Cpx1 = c("A", "B", "C", "D"),
#'   Cpx2 = c("A", "B", "C", "E"), # similar to Cpx1
#'   Cpx3 = c("F", "G", "H")
#' )
#' enrichments <- list(
#'   Cpx1 = data.frame(ID = c("GO:1", "GO:2")),
#'   Cpx2 = data.frame(ID = c("GO:1", "GO:3")), # functionally similar to Cpx1
#'   Cpx3 = data.frame(ID = c("GO:4"))
#' )
#'
#' # --- Build Network (using 2 cores for the example) ---
#' network <- buildComplexNetwork(
#'   complexes, enrichments, mode = "combined", nCores = 2
#' )
#' print(network)
#'
#' @export
#'
buildComplexNetwork <- function(complexes, enrichments,
                                mode="combined",
                                similarityMethod="jaccard",
                                alpha=0.5,
                                nCores=NULL,
                                chunkSize=1000,
                                verbose=TRUE) {
  
  if (verbose) {
    message(
      sprintf("Building complex network using '%s' similarity...",
              similarityMethod)
    )
  }
  
  # --- BUGFIX START ---
  # Handle edge case where there are not enough complexes to form pairs
  if (length(complexes) < 2) {
    if (verbose) {
      message("Fewer than 2 complexes remain; skipping network construction.")
    }
    # Return an empty tibble with the correct final column structure
    return(tibble::tibble(
      source_complex_id = character(),
      target_complex_id = character(),
      compSim = double(),
      funcSim = double(),
      sharedProt = integer(),
      sharedFunc = integer(),
      weight = double(),
      similarity_mode = character()
    ))
  }
  # --- BUGFIX END ---
  
  # --- 1. Setup Parallel Processing ---
  if (is.null(nCores)) {
    nCores <- max(1, parallel::detectCores() - 1)
  }
  if (verbose) {
    message(sprintf("Using %d cores for parallel processing.", nCores))
  }
  
  oldPlan <- future::plan(future::multisession, workers=nCores)
  on.exit(future::plan(oldPlan), add=TRUE)
  
  # --- 2. Prepare Data for Pairwise Comparison ---
  complexIds <- names(complexes)
  nComplexes <- length(complexIds)
  nPairs <- choose(nComplexes, 2)
  
  if (verbose) message(sprintf("Processing %d complex pairs...", nPairs))
  
  pairIndices <- utils::combn(nComplexes, 2, simplify=FALSE)
  
  # ... (rest of the function remains identical)
  
  # Split pairs into chunks for efficient parallel processing
  if (length(pairIndices) > 0) {
    # Create a grouping factor for splitting
    nChunks <- ceiling(length(pairIndices) / chunkSize)
    groupFactor <- rep(
      seq_len(nChunks), each=chunkSize, length.out=length(pairIndices)
    )
    pairChunks <- split(pairIndices, groupFactor)
  } else {
    pairChunks <- list()
  }
  
  # --- 3. Define the Core Processing Function for One Chunk ---
  processChunk <- function(chunkPairs) {
    chunkResults <- vector("list", length(chunkPairs))
    
    for (idx in seq_along(chunkPairs)) {
      pairIdx <- chunkPairs[[idx]]
      i <- pairIdx[1]
      j <- pairIdx[2]
      c1Id <- complexIds[i]
      c2Id <- complexIds[j]
      
      # Calculate compositional and functional similarity
      compSim <- .calculateSimilarity(
        complexes[[c1Id]], complexes[[c2Id]], method=similarityMethod
      )
      funcSim <- .calculateSimilarity(
        enrichments[[c1Id]]$ID, enrichments[[c2Id]]$ID,
        method=similarityMethod
      )
      
      # Store as a simple list for memory efficiency
      sharedProt <- length(intersect(complexes[[c1Id]],
                                     complexes[[c2Id]]))
      sharedFunc <- length(intersect(enrichments[[c1Id]]$ID,
                                     enrichments[[c2Id]]$ID))
      
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
  
  # --- 4. Run Parallel Computation ---
  if (verbose && length(pairChunks) > 0) {
    message(
      sprintf("Split into %d chunks of up to %d pairs each.",
              nChunks, chunkSize)
    )
  }
  
  # Use progressr if available, otherwise run silently
  if (verbose && requireNamespace("progressr", quietly=TRUE)) {
    progressr::with_progress({
      p <- progressr::progressor(steps=nChunks)
      chunkList <- future.apply::future_lapply(pairChunks, function(chunk) {
        result <- processChunk(chunk)
        p() # Update progress
        result
      }, future.seed=TRUE)
    })
  } else {
    chunkList <- future.apply::future_lapply(
      pairChunks, processChunk, future.seed=TRUE
    )
  }
  
  # --- 5. Combine and Finalize Results ---
  if (verbose) message("Combining results from chunks...")
  
  if (length(chunkList) > 0) {
    # Flatten the list of lists and bind into a single data frame
    allResults <- unlist(chunkList, recursive=FALSE)
    finalNetworkDf <- dplyr::bind_rows(allResults)
  } else {
    # Create an empty tibble with correct columns if no pairs were processed
    finalNetworkDf <- tibble::tibble(
      source_complex_id = character(),
      target_complex_id = character(),
      compSim = double(),
      funcSim = double(),
      sharedProt = integer(),
      sharedFunc = integer()
    )
  }
  
  if (verbose) message("Calculating final weights and filtering...")
  
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
    message(
      sprintf("Network construction complete: %d edges retained.",
              nrow(finalNetworkDf))
    )
  }
  
  return(finalNetworkDf)
}