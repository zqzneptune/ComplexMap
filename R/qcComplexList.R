#' Perform Quality Control on a List of Protein Complexes
#'
#' @description
#' A diagnostic tool to assess complex size distribution and distinguish between
#' true redundancy (synonyms) and biological hierarchy (subsets).
#'
#' @details
#' **Systems Biology Rationale:**
#' To maintain a diverse functional landscape, it is critical to distinguish
#' between two types of similarity:
#' 
#' 1.  **Jaccard Similarity (Redundancy):** Complex A and B are nearly identical.
#'     *   *Action:* Merge these (using `refineComplexList`).
#' 2.  **Overlap/Simpson Similarity (Hierarchy):** Complex A is a subset of Complex B.
#'     *   *Action:* **Keep these.** Do not merge. These often represent distinct
#'         biological states (e.g., Core Complex vs. Holo-Complex).
#' 
#' This report calculates both metrics to guide your refinement strategy.
#'
#' @param complexList A list of character vectors (protein complexes).
#' @param redundancyThreshold Threshold for Jaccard similarity warning. Defaults to 0.9.
#' @param subsetThreshold Threshold for Simpson coefficient (subset) reporting. Defaults to 0.9.
#' @param verbose Logical.
#'
#' @return Invisibly returns the input `complexList`.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
qcComplexList <- function(complexList, 
                          redundancyThreshold = 0.9,
                          subsetThreshold = 0.9,
                          verbose = TRUE) {
  
  if (!verbose) return(invisible(complexList))
  
  message("\n--- Quality Control & Diversity Report ---")
  
  numComplexes <- length(complexList)
  if (numComplexes == 0) {
    message("! WARNING: Input list is empty.")
    return(invisible(complexList))
  }
  
  allProteins <- unique.default(unlist(complexList, use.names = FALSE))
  complexSizes <- lengths(complexList)
  
  # --- 1. Size Distribution ---
  message("\n[1] Compositional Complexity")
  message(sprintf("  \u2713 Complexes: %d | Proteins: %d", numComplexes, length(allProteins)))
  message(sprintf("  \u2713 Sizes: Min=%d, Med=%d, Max=%d", 
                  min(complexSizes), round(stats::median(complexSizes)), max(complexSizes)))
  
  smallCpx <- sum(complexSizes < 3)
  if (smallCpx > 0) {
    message(sprintf("  ! NOTE: %d complexes have <3 members (potentially unstable).", smallCpx))
  }
  
  # --- 2. Redundancy vs Hierarchy Analysis ---
  message("\n[2] Redundancy (Jaccard) vs. Hierarchy (Subset)")
  
  if (numComplexes < 2) {
    message("  - Skipping pairwise analysis (N < 2).")
    return(invisible(complexList))
  }
  
  # Matrix construction
  proteinIndex <- match(unlist(complexList), allProteins)
  i_indices <- rep(seq_along(complexList), complexSizes)
  
  M <- Matrix::sparseMatrix(
    i = i_indices, j = proteinIndex, x = 1,
    dims = c(numComplexes, length(allProteins))
  )
  
  Intersection <- Matrix::tcrossprod(M) # A intersect B
  
  # A. Jaccard (A int B) / (A union B)
  Union <- outer(complexSizes, complexSizes, "+") - Intersection
  Jaccard <- Intersection / Union
  diag(Jaccard) <- 0
  
  # B. Simpson (A int B) / min(A, B) -- Detects subsets
  MinSize <- outer(complexSizes, complexSizes, pmin)
  Simpson <- Intersection / MinSize
  diag(Simpson) <- 0
  
  # Analysis
  j_scores <- Jaccard[upper.tri(Jaccard)]
  s_scores <- Simpson[upper.tri(Simpson)]
  
  n_redundant <- sum(j_scores >= redundancyThreshold)
  n_subsets <- sum(s_scores >= subsetThreshold & j_scores < redundancyThreshold)
  
  # Report Redundancy
  if (n_redundant > 0) {
    message(sprintf(
      "  ! WARNING: Found %d pairs with Jaccard >= %.2f (True Redundancy).",
      n_redundant, redundancyThreshold
    ))
    message("    -> Recommendation: Merge these using refineComplexList().")
  } else {
    message("  \u2713 No high redundancy detected.")
  }
  
  # Report Hierarchy
  if (n_subsets > 0) {
    message(sprintf(
      "  i INFO: Found %d pairs with Subset Overlap >= %.2f but Low Jaccard.",
      n_subsets, subsetThreshold
    ))
    message("    -> Recommendation: KEEP distinct. These are likely biological variants (Subsets/Supersets).")
    message("    -> Merging these would decrease landscape diversity.")
  }
  
  message("\n--- QC Complete ---\n")
  return(invisible(complexList))
}