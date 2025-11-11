#' Perform Quality Control on a List of Protein Complexes
#'
#' @description
#' This function performs a quality control analysis on a list of protein
#' complexes. It provides basic statistics, analyzes the distribution of complex
#' sizes, and calculates pairwise redundancy using the Jaccard similarity index.
#'
#' @details
#' The QC process involves three main steps:
#' 1.  **Basic Statistics:** Reports the total number of complexes and unique
#'     proteins in the list.
#' 2.  **Size Distribution:** Provides a summary of the number of proteins per
#'     complex and warns if complexes have fewer than three members.
#' 3.  **Redundancy Analysis:** A sparse binary membership matrix is constructed
#'     to efficiently calculate the Jaccard similarity for all unique pairs of
#'     complexes. A warning is issued if a significant portion of pairs
#'     exceeds the specified `redundancyThreshold`.
#'
#' @param complexList A list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param redundancyThreshold A numeric value between 0 and 1. A warning will be
#'   issued for any pair of complexes with a Jaccard similarity score greater
#'   than or equal to this threshold. Defaults to 0.8.
#' @param verbose A logical value indicating whether to print progress messages
#'   and summaries to the console. Defaults to `TRUE`.
#'
#' @return
#' Invisibly returns the original `complexList` object, allowing it to be used
#' in a pipeline.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @seealso
#' `refine_complex_list()` for a method to merge redundant complexes.
#'
#' @examples
#' # Create a sample list of protein complexes
#' complex1 <- c("A", "B", "C", "D")
#' complex2 <- c("A", "B", "C", "E") # Highly redundant with complex1
#' complex3 <- c("F", "G", "H")
#' complex4 <- c("I", "J")          # Small complex
#' complex5 <- c("X", "Y", "Z")
#' sampleList <- list(
#'   C1 = complex1, C2 = complex2, C3 = complex3, C4 = complex4, C5 = complex5
#' )
#'
#' # Run the quality control analysis
#' qcComplexList(sampleList)
#'
#' @export
#'
qcComplexList <- function(complexList, redundancyThreshold=0.8,
                          verbose=TRUE) {
  
  if (verbose) {
    message("\n--- Running Quality Control on Complex List ---")
  }
  
  numComplexes <- length(complexList)
  
  if (numComplexes == 0) {
    warning("Input complex list is empty.")
    return(invisible(complexList))
  }
  
  allProteins <- unique.default(unlist(complexList, use.names=FALSE))
  numProteins <- length(allProteins)
  complexSizes <- lengths(complexList)
  
  if (verbose) {
    message("\n[1] Basic Statistics:")
    message(sprintf("    - Total number of complexes: %d", numComplexes))
    message(sprintf("    - Total number of unique proteins: %d",
                    numProteins))
    message("\n[2] Complex Size Distribution:")
    summaryStats <- utils::capture.output(summary(complexSizes))
    for (line in summaryStats) {
      message("    ", line)
    }
  }
  
  smallComplexes <- sum(complexSizes < 3)
  if (smallComplexes > 0) {
    warning(
      sprintf("%d complexes have fewer than 3 members.", smallComplexes)
    )
  }
  
  if (verbose) message("\n[3] Redundancy Analysis:")
  
  if (numComplexes < 2) {
    if (verbose) message("    - Skipping: not enough complexes to analyze.")
  } else {
    proteinIndex <- stats::setNames(seq_along(allProteins), allProteins)
    i_indices <- rep(seq_along(complexList), times=complexSizes)
    j_indices <- proteinIndex[unlist(complexList, use.names=FALSE)]
    
    membershipMatrix <- Matrix::sparseMatrix(
      i=i_indices, j=j_indices, x=1,
      dims=c(numComplexes, numProteins)
    )
    
    intersectionMatrix <- membershipMatrix %*% Matrix::t(membershipMatrix)
    sizeSumMatrix <- outer(complexSizes, complexSizes, "+")
    unionMatrix <- sizeSumMatrix - intersectionMatrix
    
    jaccardMatrix <- intersectionMatrix / unionMatrix
    diag(jaccardMatrix) <- 0
    simScores <- jaccardMatrix[upper.tri(jaccardMatrix)]
    
    if (verbose) {
      message("    - Distribution of Jaccard similarity scores:")
      summaryStats <- utils::capture.output(summary(simScores))
      for (line in summaryStats) {
        message("    ", line)
      }
    }
    
    numRedundant <- sum(simScores >= redundancyThreshold, na.rm=TRUE)
    if (numRedundant > 0) {
      nPairs <- length(simScores)
      percentRedundant <- (numRedundant / nPairs) * 100
      warning(
        sprintf(
          "%d pairs (%.2f%%) are highly redundant (Jaccard >= %.2f).",
          numRedundant, percentRedundant, redundancyThreshold
        ),
        "\n    Consider using a refinement function to merge them."
      )
    } else if (verbose) {
      message("    - No highly redundant complex pairs detected.")
    }
  }
  
  if (verbose) message("\n--- QC Complete ---\n")
  return(invisible(complexList))
}