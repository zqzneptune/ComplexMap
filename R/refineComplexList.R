#' Refine a List of Protein Complexes by Size and Redundancy
#'
#' @description
#' This function refines a list of protein complexes by first filtering them
#' based on size and then merging highly redundant complexes based on a Jaccard
#' similarity threshold.
#'
#' @details
#' The refinement process consists of two main stages:
#' 1.  **Size Filtering:** Complexes that are smaller than `minSize` or larger
#'     than `maxSize` are removed from the list.
#' 2.  **Redundancy Merging:** A Jaccard similarity matrix is calculated for all
#'     remaining pairs of complexes. A union-find algorithm is then used to
#'     identify clusters (or "merge groups") of complexes where all members are
#'     connected by a similarity score greater than or equal to the
#'     `mergeThreshold`. These groups are then merged into single complexes.
#'
#' Finally, all complexes in the refined list are renamed to a standardized
#' format ("CpxMap_0001", "CpxMap_0002", etc.).
#'
#' @param complexList A named list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param minSize An integer specifying the minimum number of proteins a complex
#'   must have to be retained. Defaults to 3.
#' @param maxSize An integer specifying the maximum number of proteins a complex
#'   can have to be retained. Defaults to 500.
#' @param mergeThreshold A numeric value (0-1) for the Jaccard similarity.
#'   Complexes with a score >= this value will be merged. Defaults to 0.9.
#' @param verbose A logical value indicating whether to print progress messages.
#'   Defaults to `TRUE`.
#'
#' @return
#' A refined and renamed named list of protein complexes.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @seealso
#' `qcComplexList()` for quality control analysis of a complex list.
#'
#' @examples
#' # Create a sample list of protein complexes
#' complex1 <- c("A", "B", "C", "D")
#' complex2 <- c("A", "B", "C", "E") # Highly redundant with complex1
#' complex3 <- c("F", "G", "H")
#' complex4 <- c("I", "J")          # Too small, will be filtered
#' sampleList <- list(
#'   C1 = complex1, C2 = complex2, C3 = complex3, C4 = complex4
#' )
#'
#' # Refine the list using a high similarity threshold
#' refinedList <- refineComplexList(sampleList, mergeThreshold = 0.75)
#'
#' @export
#'
refineComplexList <- function(complexList, minSize=3, maxSize=500,
                              mergeThreshold=0.9, verbose=TRUE) {
  
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
    return(complexListFiltered)
  }
  
  # 2. Merge using original union-find logic
  if (verbose) {
    message(sprintf("Identifying merge groups with Jaccard >= %.2f...",
                    mergeThreshold))
  }
  
  nComplexes <- length(complexListFiltered)
  complexIds <- names(complexListFiltered)
  
  # Efficiently build Jaccard similarity matrix
  allProteins <- unique.default(unlist(complexListFiltered, use.names=FALSE))
  proteinIndex <- stats::setNames(seq_along(allProteins), allProteins)
  complexSizesFiltered <- lengths(complexListFiltered)
  
  i_indices <- rep(seq_along(complexListFiltered), times=complexSizesFiltered)
  j_indices <- proteinIndex[unlist(complexListFiltered, use.names=FALSE)]
  
  membershipMatrix <- Matrix::sparseMatrix(
    i=i_indices, j=j_indices, x=1,
    dims=c(nComplexes, length(allProteins))
  )
  intersectionMatrix <- Matrix::tcrossprod(membershipMatrix)
  sizeSumMatrix <- outer(complexSizesFiltered, complexSizesFiltered, "+")
  unionMatrix <- sizeSumMatrix - intersectionMatrix
  simMatrix <- intersectionMatrix / unionMatrix
  diag(simMatrix) <- 0
  
  simMatrixDense <- as.matrix(simMatrix)
  upperTriIndices <- which(
    upper.tri(simMatrixDense) & simMatrixDense >= mergeThreshold,
    arr.ind=TRUE
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
  
  resultList <- vector("list", nMerges)
  clusterGroups <- split(seq_len(nComplexes), clusters)
  
  idx <- 1
  for (clusterMembers in clusterGroups) {
    memberIds <- complexIds[clusterMembers]
    if (length(clusterMembers) == 1) {
      resultList[[idx]] <- complexListFiltered[[clusterMembers[1]]]
      names(resultList)[idx] <- memberIds
    } else {
      mergedProteins <- unique.default(
        unlist(complexListFiltered[clusterMembers], use.names=FALSE)
      )
      if (length(memberIds) <= 3) {
        newId <- paste(memberIds, collapse="_")
      } else {
        newId <- sprintf("%s_%s_plus%d_merged", memberIds[1],
                         memberIds[2], length(memberIds) - 2)
      }
      resultList[[idx]] <- mergedProteins
      names(resultList)[idx] <- newId
    }
    idx <- idx + 1
  }
  
  if (verbose) {
    message(sprintf("Merging complete. Final list has %d complexes.",
                    length(resultList)))
  }
  
  names(resultList) <- sprintf("CpxMap_%04d", seq_along(resultList))
  
  if (verbose) message("\n--- Refinement Complete ---\n")
  
  return(resultList)
}