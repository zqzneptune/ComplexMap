#' Perform Quality Control on a List of Protein Complexes
#'
#' @description
#' This function performs a quality control analysis on a list of protein
#' complexes, delivering a user-friendly report with key statistics and
#' actionable warnings.
#'
#' @details
#' The QC process involves three main steps, presented in a clear report:
#' 1.  **Basic Statistics:** Reports the total number of complexes and unique
#'     proteins.
#'     
#' 2.  **Size Distribution:** Summarizes the number of proteins per complex,
#'     highlighting the smallest, median, and largest complexes. It will
#'     issue an in-report warning if complexes have fewer than 3 members.
#'     
#' 3.  **Redundancy Analysis:** Calculates the Jaccard similarity for all pairs
#'     of complexes. It reports the median and maximum similarity and issues an
#'     in-report warning if any pairs exceed the `redundancyThreshold`.
#'
#' This function is intended as a **diagnostic tool**. It uses the Jaccard index
#' for its clear and intuitive interpretation of redundancy. For actively merging
#' complexes, see the more advanced `refineComplexList()` function.
#'
#' @param complexList A list where each element is a character vector of
#'   protein identifiers representing a complex.
#' @param redundancyThreshold A numeric value between 0 and 1. A warning will be
#'   issued for any pair of complexes with a Jaccard similarity score greater
#'   than or equal to this threshold. Defaults to 0.8.
#' @param verbose A logical value indicating whether to print the QC report to
#'   the console. Defaults to `TRUE`.
#'
#' @return
#' Invisibly returns the original `complexList` object, allowing it to be used
#' in a pipeline.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @seealso
#' `refineComplexList()` for a function to actively merge redundant complexes.
#' 
#' @export
#' 
#' @examples
#' # Create a sample list of protein complexes
#' complex1 <- c("A", "B", "C", "D")
#' complex2 <- c("A", "B", "C", "E") # Highly redundant with complex1
#' complex3 <- c("F", "G", "H")
#' complex4 <- c("I", "J")          # Small complex
#' sampleList <- list(
#'   C1 = complex1, C2 = complex2, C3 = complex3, C4 = complex4
#' )
#'
#' # Run the quality control analysis to see the new report format
#' qcComplexList(sampleList)
#'
qcComplexList <- function(complexList, redundancyThreshold = 0.8,
                          verbose = TRUE) {
  
  if (!verbose) {
    return(invisible(complexList))
  }
  
  message("\n--- Quality Control Report for Complex List ---")
  
  numComplexes <- length(complexList)
  
  if (numComplexes == 0) {
    message("! WARNING: Input complex list is empty. No analysis to perform.")
    message("\n--- QC Complete ---\n")
    return(invisible(complexList))
  }
  
  allProteins <- unique.default(unlist(complexList, use.names = FALSE))
  numProteins <- length(allProteins)
  complexSizes <- lengths(complexList)
  
  # --- Section 1: Basic Statistics ---
  message("\n[1] Basic Statistics")
  message(sprintf("  \u2713 Total Complexes: %d", numComplexes))
  message(sprintf("  \u2713 Unique Proteins: %d", numProteins))
  
  # --- Section 2: Complex Size Distribution ---
  message("\n[2] Complex Size Distribution")
  message(sprintf("  - Smallest Complex: %d proteins", min(complexSizes)))
  message(sprintf("  - Median Complex:   %d proteins", round(stats::median(complexSizes))))
  message(sprintf("  - Largest Complex:  %d proteins", max(complexSizes)))
  
  smallComplexesCount <- sum(complexSizes < 3)
  if (smallComplexesCount > 0) {
    message(sprintf(
      "  ! WARNING: Found %d complex(es) with fewer than 3 members.",
      smallComplexesCount
    ))
    message("    -> These are often considered too small for robust analysis.")
  } else {
    message("  \u2713 All complexes have 3 or more members.")
  }
  
  # --- Section 3: Redundancy Analysis ---
  message("\n[3] Redundancy Analysis (Jaccard Index)")
  if (numComplexes < 2) {
    message("  - Skipping: Not enough complexes to compare.")
  } else {
    proteinIndex <- stats::setNames(seq_along(allProteins), allProteins)
    i_indices <- rep(seq_along(complexList), times = complexSizes)
    j_indices <- proteinIndex[unlist(complexList, use.names = FALSE)]
    
    membershipMatrix <- Matrix::sparseMatrix(
      i = i_indices, j = j_indices, x = 1,
      dims = c(numComplexes, numProteins)
    )
    
    intersectionMatrix <- membershipMatrix %*% Matrix::t(membershipMatrix)
    sizeSumMatrix <- outer(complexSizes, complexSizes, "+")
    unionMatrix <- sizeSumMatrix - intersectionMatrix
    
    jaccardMatrix <- intersectionMatrix / unionMatrix
    diag(jaccardMatrix) <- 0
    simScores <- jaccardMatrix[upper.tri(jaccardMatrix)]
    
    message(sprintf("  - Median Similarity: %.3f", stats::median(simScores, na.rm = TRUE)))
    message(sprintf("  - Max Similarity:    %.3f", max(simScores, na.rm = TRUE)))
    
    numRedundant <- sum(simScores >= redundancyThreshold, na.rm = TRUE)
    if (numRedundant > 0) {
      nPairs <- length(simScores)
      percentRedundant <- (numRedundant / nPairs) * 100
      message(sprintf(
        "  ! WARNING: Found %d pairs (%.2f%%) with similarity >= %.2f.",
        numRedundant, percentRedundant, redundancyThreshold
      ))
      message("    -> Consider using refineComplexList() to merge them.")
    } else {
      message(sprintf(
        "  \u2713 No highly redundant pairs detected (>= %.2f).",
        redundancyThreshold
      ))
    }
  }
  
  message("\n--- QC Complete ---\n")
  return(invisible(complexList))
}