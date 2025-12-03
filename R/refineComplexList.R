#' Refine a List of Protein Complexes by Size and Redundancy
#'
#' @description
#' This function filters complexes by size and merges highly redundant ones.
#' 
#' @details
#' **Systems Biology Rationale:**
#' To preserve the functional diversity of the input landscape, this function
#' defaults to a conservative **Jaccard** similarity with a high threshold (0.90).
#' 
#' * **Trust the Input:** We assume the upstream clustering method identified
#'   biologically relevant variants (e.g., a complex with vs. without a regulatory subunit).
#' * **Minimal Merging:** We only merge complexes if they are effectively
#'   synonyms (nearly identical composition). 
#'   
#' We discourage using "overlap" or "matching_score" for this step, as they 
#' tend to merge sub-complexes into parents, destroying the specific functional 
#' signals you wish to visualize.
#'
#' @param complexList A named list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param minSize Minimum complex size. Defaults to 3.
#' @param maxSize Maximum complex size. Defaults to 500.
#' @param mergeThreshold Numeric (0-1). Threshold for merging. Defaults to **0.90** 
#'   (strict) to preserve distinct variants.
#' @param similarityMethod Metric for merging. Defaults to **"jaccard"**.
#'   Options:
#'   * `"jaccard"`: (Default) `Intersection / Union`. Penalizes size differences. 
#'     Best for preserving specific variants.
#'   * `"overlap"`: `Intersection / Min(A,B)`. Merges subsets into parents. 
#'     Use only if you want to aggressively reduce data.
#'   * `"matching_score"`: `Intersection^2 / (SizeA * SizeB)`. Geometric approach.
#'   * `"dice"`: `2 * Intersection / (SizeA + SizeB)`.
#' @param verbose Logical.
#'
#' @return
#' A list containing:
#' * `refinedComplexes`: The final list of merged complexes.
#' * `mergeMap`: A tibble mapping original IDs to final IDs.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
refineComplexList <- function(complexList, minSize = 3, maxSize = 500,
                              mergeThreshold = 0.90,
                              similarityMethod = "jaccard", verbose = TRUE) {
  
  if (verbose) message("\n--- Refining Input Complex List (Minimal Merging Strategy) ---")
  
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
    if (verbose) message("Not enough complexes to merge. Returning size-filtered list.")
    finalIds <- sprintf("CpxMap_%04d", seq_along(complexListFiltered))
    names(complexListFiltered) <- finalIds
    mergeMap <- tibble::tibble(originalId = names(complexListFiltered), finalId = finalIds)
    return(list(refinedComplexes = complexListFiltered, mergeMap = mergeMap))
  }
  
  # 2. Merge using union-find logic
  if (verbose) {
    message(sprintf("Identifying merge groups: method='%s', threshold >= %.2f...",
                    similarityMethod, mergeThreshold))
  }
  
  nComplexes <- length(complexListFiltered)
  originalComplexIds <- names(complexListFiltered)
  sizes <- lengths(complexListFiltered)
  
  # Efficiently build the intersection matrix
  allProteins <- unique.default(unlist(complexListFiltered, use.names = FALSE))
  proteinIndex <- match(unlist(complexListFiltered, use.names = FALSE), allProteins)
  i_indices <- rep(seq_along(complexListFiltered), times = sizes)
  
  membershipMatrix <- Matrix::sparseMatrix(
    i = i_indices, j = proteinIndex, x = 1,
    dims = c(nComplexes, length(allProteins))
  )
  
  intersectionMatrix <- Matrix::tcrossprod(membershipMatrix)
  
  simMatrix <- switch(
    similarityMethod,
    "jaccard" = {
      sizeSumMatrix <- outer(sizes, sizes, "+")
      unionMatrix <- sizeSumMatrix - intersectionMatrix
      intersectionMatrix / unionMatrix
    },
    "overlap" = {
      minSizeMatrix <- outer(sizes, sizes, pmin)
      intersectionMatrix / minSizeMatrix
    },
    "matching_score" = {
      prodSizeMatrix <- outer(sizes, sizes, "*")
      (intersectionMatrix^2) / prodSizeMatrix
    },
    "dice" = {
      sumSizeMatrix <- outer(sizes, sizes, "+")
      (2 * intersectionMatrix) / sumSizeMatrix
    },
    stop("Invalid 'similarityMethod'. Choose 'jaccard', 'overlap', 'matching_score', or 'dice'.")
  )
  
  simMatrix[is.nan(simMatrix)] <- 0
  diag(simMatrix) <- 0
  
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
  
  if (verbose) {
    nMergedGroups <- sum(tabulate(clusters) > 1)
    message(sprintf("Found %d redundancy groups. Merging %d complexes into %d.",
                    nMergedGroups, nComplexes, nMerges))
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
  if (verbose) message("\n--- Refinement Complete ---\n")
  return(list(refinedComplexes = refinedComplexes, mergeMap = mergeMap))
}