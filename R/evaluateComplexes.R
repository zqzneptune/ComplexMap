#' Evaluate Predicted Protein Complexes Against a Reference Set
#'
#' @description
#' Calculates four standard metrics for evaluating protein complex predictions:
#' Positive Predictive Value (PPV), Sensitivity (Sn), Accuracy (Acc), and the
#' Maximum Matching Ratio (MMR).
#'
#' @details
#' This function is optimized for speed by calculating a shared intersection
#' matrix between predicted and reference complexes in parallel. This matrix is
#' then used as the basis for all four metric calculations.
#'
#' - **PPV, Sn, and Acc** are calculated based on the confusion matrix between
#'   predicted and reference complexes, as described in the literature (e.g.,
#'   Zhang et al., 2012).
#' 
#' - **MMR** is calculated by first deriving an overlap score matrix, where the
#'   score for a predicted complex (P) and a reference complex (R) is
#'   `|P ∩ R|² / (|P| * |R|)`. The
#'   [Hungarian algorithm](https://en.wikipedia.org/wiki/Hungarian_algorithm)
#'   is then used to solve the maximum weight bipartite matching problem. This
#'   approach is based on the method described by Nepusz et al. (2012).
#'
#' The parallel computation uses `parallel::mclapply`, which is not available on
#' Windows. On Windows, the calculation will run sequentially.
#'
#' @references
#' Nepusz, T., Yu, H. & Paccanaro, A. (2012). Detecting overlapping protein
#' complexes in protein-protein interaction networks. *Nature Methods*, 9,
#' 471–472. \doi{10.1038/nmeth.1938}
#'
#' Zhang XF, Dai DQ, Ou-Yang L, Wu MY (2012). Exploring Overlapping Functional
#' Units with Various Structure in Protein Interaction Networks. *PLOS ONE*,
#' 7(8): e43092. \doi{10.1371/journal.pone.0043092}
#'
#' @param predictedComplexes A list of predicted protein complexes.
#' @param referenceComplexes A list of reference (gold standard) complexes.
#' @param nCores The number of CPU cores to use for parallel computation.
#'   Defaults to one less than the total number of detected cores.
#' @param verbose A logical value indicating whether to print progress messages.
#'   Defaults to `TRUE`.
#'
#' @return
#' A named list containing four numeric values: `PPV`, `Sn`, `Acc`, and `MMR`.
#' Returns `NA` for all metrics if either input list is empty.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data ---
#' # Predicted complexes
#' pred1 <- c("A", "B", "C")
#' pred2 <- c("D", "E", "F")
#' pred3 <- c("A", "G", "H")
#' predicted <- list(P1 = pred1, P2 = pred2, P3 = pred3)
#'
#' # Reference complexes (gold standard)
#' ref1 <- c("A", "B", "C", "X") # Good match for pred1
#' ref2 <- c("D", "E", "F")     # Perfect match for pred2
#' ref3 <- c("I", "J", "K")     # Unmatched complex
#' reference <- list(R1 = ref1, R2 = ref2, R3 = ref3)
#'
#' # --- Evaluation ---
#' # Use 2 cores for the example
#' metrics <- evaluateComplexes(predicted, reference, nCores = 2)
#' print(metrics)
#'
#' @export
evaluateComplexes <- function(predictedComplexes, referenceComplexes,
                              nCores=NULL, verbose=TRUE) {
  if (verbose) message("\n--- Evaluating Complex Predictions ---")
  
  if (is.null(nCores)) {
    nCores <- max(1, parallel::detectCores() - 1)
  }
  
  numPredicted <- length(predictedComplexes)
  numReference <- length(referenceComplexes)
  
  if (numReference == 0 || numPredicted == 0) {
    warning("One or both complex lists are empty.")
    return(list(PPV=NA, Sn=NA, Acc=NA, MMR=NA))
  }
  
  # Check for any overlap between protein sets
  allRefProteins <- unique.default(unlist(referenceComplexes, use.names=FALSE))
  allPredProteins <- unique.default(unlist(predictedComplexes, use.names=FALSE))
  
  if (length(intersect(allRefProteins, allPredProteins)) == 0) {
    if (verbose) message("No overlap detected between protein sets.")
    return(list(PPV=0, Sn=0, Acc=0, MMR=0))
  }
  
  # Map each protein to its reference complexes
  refProteinMap <- split(
    rep(seq_along(referenceComplexes), lengths(referenceComplexes)),
    unlist(referenceComplexes, use.names=FALSE)
  )
  
  # --- SHARED COMPUTATION: Calculate intersection counts in parallel ---
  if (verbose) {
    message(
      sprintf("[1] Calculating intersections using %d core(s)...", nCores)
    )
  }
  intersectionList <- parallel::mclapply(predictedComplexes, function(predSet) {
    candidateRefs <- unique.default(unlist(
      refProteinMap[intersect(predSet, names(refProteinMap))],
      use.names=FALSE
    ))
    
    intersections <- integer(numReference)
    if (length(candidateRefs) > 0) {
      for (j in candidateRefs) {
        intersections[j] <- length(
          intersect(predSet, referenceComplexes[[j]])
        )
      }
    }
    intersections
  }, mc.cores=nCores)
  
  # Rows = predicted, Cols = reference
  intersectionMatrix <- do.call(rbind, intersectionList)
  
  # --- METRICS 1, 2, 3: Calculate Sn, PPV, Acc ---
  if (verbose) message("[2] Calculating PPV, Sensitivity, and Accuracy...")
  refSizes <- lengths(referenceComplexes)
  predSizesInRef <- lengths(lapply(predictedComplexes,
                                   intersect, allRefProteins))
  
  maxPerRef <- apply(intersectionMatrix, 2, max)
  Sn <- sum(maxPerRef) / sum(refSizes)
  
  maxPerPred <- apply(intersectionMatrix, 1, max)
  totalSizePredInRef <- sum(predSizesInRef)
  PPV <- if (totalSizePredInRef > 0) sum(maxPerPred)/totalSizePredInRef else 0
  
  Acc <- sqrt(abs(Sn * PPV))
  
  # --- METRIC 4: Calculate MMR ---
  if (verbose) message("[3] Calculating Maximum Matching Ratio (MMR)...")
  predSizes <- lengths(predictedComplexes)
  
  # Vectorized calculation of the overlap matrix
  overlapMatrix <- (intersectionMatrix^2) / outer(predSizes, refSizes)
  overlapMatrix[is.nan(overlapMatrix)] <- 0
  
  # Pad the smaller dimension with zeros for the solver
  maxDim <- max(numPredicted, numReference)
  if (nrow(overlapMatrix) != maxDim || ncol(overlapMatrix) != maxDim) {
    paddedMatrix <- matrix(0, nrow=maxDim, ncol=maxDim)
    paddedMatrix[seq_len(numPredicted), seq_len(numReference)] <- overlapMatrix
    overlapMatrix <- paddedMatrix
  }
  
  # Solve maximum weighted bipartite matching
  assignment <- clue::solve_LSAP(overlapMatrix, maximum=TRUE)
  
  # Vectorized calculation of total weight for valid assignments
  validRows <- seq_len(numPredicted)[assignment <= numReference]
  matchedPairs <- cbind(validRows, assignment[validRows])
  totalWeight <- sum(overlapMatrix[matchedPairs])
  
  MMR <- totalWeight / numReference
  
  if (verbose) message("\n--- Evaluation Complete ---\n")
  
  return(list(PPV=PPV, Sn=Sn, Acc=Acc, MMR=MMR))
}