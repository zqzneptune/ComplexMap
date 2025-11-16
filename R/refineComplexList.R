#' Refine a List of Protein Complexes by Size and Redundancy
#'
#' @description
#' This function refines a list of protein complexes by first filtering them
#' based on size and then merging highly redundant complexes based on a
#' similarity threshold.
#'
#' @details
#' The refinement process consists of two main stages:
#' 
#' 1.  **Size Filtering:** Complexes smaller than `minSize` or larger than
#'     `maxSize` are removed.
#'     
#' 2.  **Redundancy Merging:** A similarity matrix is calculated for all
#'     remaining complex pairs using the chosen `similarityMethod`. A
#'     union-find algorithm identifies clusters of complexes connected by a
#'     similarity score >= `mergeThreshold`, which are then merged.
#'
#' Finally, all complexes in the refined list are renamed to a standardized
#' format ("CpxMap_0001", etc.), and a traceability map is generated.
#'
#' @param complexList A named list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param minSize An integer specifying the minimum number of proteins a complex
#'   must have to be retained. Defaults to 3.
#' @param maxSize An integer specifying the maximum number of proteins a complex
#'   can have to be retained. Defaults to 500.
#' @param mergeThreshold A numeric value (0-1). Complexes with a similarity
#'   score >= this value will be merged. Defaults to 0.9.
#' @param similarityMethod A character string specifying the similarity metric
#'   for merging. Defaults to `"matching_score"`. The available options are:
#'   \describe{
#'     \item{`"matching_score"`}{`Intersection² / (|A| * |B|)`. Rewards large,
#'       shared cores. Aligns with the MMR evaluation metric.}
#'     \item{`"simpson"`}{`Intersection / min(|A|, |B|)`. Also known as the
#'       Overlap coefficient. excels at merging sub-complexes into larger parents.}
#'     \item{`"jaccard"`}{`Intersection / Union`. A classic, balanced metric
#'       that penalizes size differences.}
#'     \item{`"dice"`}{`2 * Intersection / (|A| + |B|)`. Similar to Jaccard but
#'       generally less stringent.}
#'   }
#' @param verbose A logical value indicating whether to print progress messages.
#'   Defaults to `TRUE`.
#'
#' @return
#' A list containing two named elements:
#' \describe{
#'   \item{`refinedComplexes`}{The final, renamed list of merged complexes.}
#'   \item{`mergeMap`}{A `tibble` with columns `originalId` and `finalId`,
#'   mapping each original complex to its final standardized ID.}
#' }
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Create a sample list of protein complexes
#' c1 <- c("A", "B", "C", "D", "E")
#' c2 <- c("A", "B", "C", "D", "F") # High similarity with c1
#' c3 <- c("A", "B", "C")          # Subset of c1 and c2
#' c4 <- c("X", "Y", "Z")
#' sampleList <- list(C1=c1, C2=c2, C3=c3, C4=c4)
#'
#' # Refine using the default "matching_score"
#' refineComplexList(sampleList, mergeThreshold = 0.6)
#'
#' # Refine using the "simpson" method to merge the subset
#' refineComplexList(sampleList, mergeThreshold = 0.9, similarityMethod = "simpson")
#'
refineComplexList <- function(complexList, minSize = 3, maxSize = 500,
                              mergeThreshold = 0.9,
                              similarityMethod = "matching_score", verbose = TRUE) {
  
  if (verbose) message("\n--- Refining Input Complex List ---")
  
  # 1. Filter by size
  initialCount <- length(complexList)
  complexSizes <- lengths(complexList)
  keepIdx <- complexSizes >= minSize & complexSizes <= maxSize
  complexListFiltered <- complexList[keepIdx]
  nFiltered <- initialCount - length(complexListFiltered)
  
  if (verbose) {
    message(sprintf("Filtered %d complexes by size. Retaining %d.",
                    nFiltered, length(complexListFiltered)))
  }
  
  if (length(complexListFiltered) < 2) {
    if (verbose) {
      message("Not enough complexes to merge. Returning size-filtered list.")
    }
    finalIds <- sprintf("CpxMap_%04d", seq_along(complexListFiltered))
    names(complexListFiltered) <- finalIds
    mergeMap <- tibble::tibble(
      originalId = names(complexListFiltered),
      finalId = finalIds
    )
    return(list(refinedComplexes = complexListFiltered, mergeMap = mergeMap))
  }
  
  # 2. Merge using union-find logic
  if (verbose) {
    message(sprintf("Identifying merge groups with %s >= %.2f...",
                    similarityMethod, mergeThreshold))
  }
  
  nComplexes <- length(complexListFiltered)
  originalComplexIds <- names(complexListFiltered)
  sizes <- lengths(complexListFiltered)
  
  # Efficiently build the intersection matrix (shared foundation for all metrics)
  allProteins <- unique.default(unlist(complexListFiltered, use.names = FALSE))
  proteinIndex <- stats::setNames(seq_along(allProteins), allProteins)
  i_indices <- rep(seq_along(complexListFiltered), times = sizes)
  j_indices <- proteinIndex[unlist(complexListFiltered, use.names = FALSE)]
  
  membershipMatrix <- Matrix::sparseMatrix(
    i = i_indices, j = j_indices, x = 1,
    dims = c(nComplexes, length(allProteins))
  )
  intersectionMatrix <- Matrix::tcrossprod(membershipMatrix)
  
  # Calculate similarity matrix based on the chosen method
  simMatrix <- switch(
    similarityMethod,
    "matching_score" = {
      (intersectionMatrix^2) / outer(sizes, sizes, "*")
    },
    "simpson" = {
      intersectionMatrix / outer(sizes, sizes, pmin)
    },
    "jaccard" = {
      sizeSumMatrix <- outer(sizes, sizes, "+")
      unionMatrix <- sizeSumMatrix - intersectionMatrix
      intersectionMatrix / unionMatrix
    },
    "dice" = {
      (2 * intersectionMatrix) / outer(sizes, sizes, "+")
    },
    stop("Invalid 'similarityMethod'. Choose 'matching_score', 'simpson', 'jaccard', or 'dice'.")
  )
  
  simMatrix[is.nan(simMatrix)] <- 0 # Handle potential division by zero
  diag(simMatrix) <- 0
  
  # --- Continue with the merging and traceability logic (no changes needed here) ---
  simMatrixDense <- as.matrix(simMatrix)
  upperTriIndices <- which(
    upper.tri(simMatrixDense) & simMatrixDense >= mergeThreshold,
    arr.ind = TRUE
  )
  
  parent <- seq_len(nComplexes)
  findRoot <- function(i) {
    path <- integer(0)
    while (parent[i] != i) {
      path <- c(path, i)
      i <- parent[i]
    }
    if (length(path) > 0) parent[path] <<- i
    return(i)
  }
  
  if (nrow(upperTriIndices) > 0) {
    for (k in seq_len(nrow(upperTriIndices))) {
      i <- upperTriIndices[k, 1]
      j <- upperTriIndices[k, 2]
      rootI <- findRoot(i)
      rootJ <- findRoot(j)
      if (rootI != rootJ) parent[rootJ] <- rootI
    }
  }
  
  clusters <- vapply(seq_len(nComplexes), findRoot, integer(1))
  uniqueClusters <- unique(clusters)
  nMerges <- length(uniqueClusters)
  nMergeOps <- nComplexes - nMerges
  
  if (verbose) {
    message(sprintf("Found %d merge groups. Merging %d complexes into %d.",
                    sum(tabulate(clusters) > 1), nMergeOps, nMerges))
  }
  
  refinedComplexes <- vector("list", nMerges)
  finalIds <- sprintf("CpxMap_%04d", seq_along(uniqueClusters))
  names(refinedComplexes) <- finalIds
  
  mergeMapList <- vector("list", nMerges)
  clusterGroups <- split(seq_len(nComplexes), clusters)
  
  for (i in seq_along(clusterGroups)) {
    clusterMembers <- clusterGroups[[i]]
    finalId <- finalIds[i]
    
    mergedProteins <- unique.default(
      unlist(complexListFiltered[clusterMembers], use.names = FALSE)
    )
    refinedComplexes[[finalId]] <- mergedProteins
    
    originalIdsInCluster <- originalComplexIds[clusterMembers]
    mergeMapList[[i]] <- tibble::tibble(
      originalId = originalIdsInCluster,
      finalId = finalId
    )
  }
  
  mergeMap <- dplyr::bind_rows(mergeMapList)
  
  if (verbose) {
    message(sprintf("Merging complete. Final list has %d complexes.",
                    length(refinedComplexes)))
  }
  if (verbose) message("\n--- Refinement Complete ---\n")
  
  return(list(refinedComplexes = refinedComplexes, mergeMap = mergeMap))
}