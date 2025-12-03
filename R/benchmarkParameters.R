#' Benchmark Complex Refinement Parameters
#'
#' @description
#' Evaluates how different merging thresholds affect the quality AND diversity
#' of the resulting complex map.
#'
#' @details
#' **Systems Biology Rationale:**
#' When optimizing the `mergeThreshold`, there is a trade-off between removing
#' redundancy (increasing PPV) and losing specific variants (decreasing Sensitivity
#' and Diversity).
#' 
#' This function returns the standard metrics (PPV, Sn, Acc, MMR) but also adds
#' **NumComplexes**.
#' - **Collapse Warning:** If `NumComplexes` drops precipitously between two thresholds,
#'   you have likely crossed the point where distinct biological variants are being
#'   merged into generic blobs.
#' - **Recommendation:** Choose the threshold that maximizes F1/MMR while maintaining
#'   a `NumComplexes` count close to your expected biological complexity.
#'
#' @param complexList A list of predicted protein complexes to be refined.
#' @param referenceComplexes A list of reference (gold standard) complexes.
#' @param threshold_range Numeric vector of thresholds to test. 
#'   Defaults to `seq(0.5, 1.0, by = 0.05)`.
#' @param similarityMethod The metric used for merging. Defaults to **"jaccard"**.
#' @param ... Additional arguments passed to `refineComplexList` (e.g., minSize).
#'
#' @return
#' A `tibble` with columns:
#' - `mergeThreshold`: The threshold tested.
#' - `NumComplexes`: Number of complexes remaining (Diversity metric).
#' - `PPV`: Positive Predictive Value (Precision).
#' - `Sn`: Sensitivity (Recall).
#' - `F1`: Harmonic mean of PPV and Sn.
#' - `Acc`: Accuracy (Geometric mean of PPV and Sn).
#' - `MMR`: Maximum Matching Ratio.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
benchmarkParameters <- function(complexList, referenceComplexes,
                                threshold_range = seq(0.5, 1.0, by = 0.05),
                                similarityMethod = "jaccard",
                                ...) {
  
  message("--- Starting Parameter Benchmarking ---")
  message(sprintf("Method: %s (Focus: Balancing Accuracy vs. Diversity)", similarityMethod))
  
  # Initialize a list to store the results
  results_list <- vector("list", length(threshold_range))
  
  for (i in seq_along(threshold_range)) {
    current_threshold <- threshold_range[i]
    
    if (interactive()) {
      message(sprintf("Testing mergeThreshold = %.2f...", current_threshold))
    }
    
    # 1. Refine
    refinement_output <- refineComplexList(
      complexList = complexList,
      mergeThreshold = current_threshold,
      similarityMethod = similarityMethod,
      verbose = FALSE,
      ...
    )
    
    refined_complexes <- refinement_output$refinedComplexes
    num_refined <- length(refined_complexes)
    
    # 2. Evaluate
    if (num_refined == 0) {
      metrics <- tibble::tibble(
        mergeThreshold = current_threshold,
        NumComplexes = 0,
        PPV = NA_real_, Sn = NA_real_, F1 = NA_real_,
        Acc = NA_real_, MMR = NA_real_
      )
    } else {
      eval_metrics <- evaluateComplexes(
        predictedComplexes = refined_complexes,
        referenceComplexes = referenceComplexes,
        verbose = FALSE
      )
      
      # Calculate F1 Score (Harmonic Mean)
      # Safe calculation to avoid division by zero
      ppv <- eval_metrics$PPV
      sn <- eval_metrics$Sn
      f1 <- if (!is.na(ppv) && !is.na(sn) && (ppv + sn) > 0) {
        2 * (ppv * sn) / (ppv + sn)
      } else {
        0
      }
      
      metrics <- tibble::tibble(
        mergeThreshold = current_threshold,
        NumComplexes = num_refined,
        PPV = ppv,
        Sn = sn,
        F1 = f1,
        Acc = eval_metrics$Acc,
        MMR = eval_metrics$MMR
      )
    }
    results_list[[i]] <- metrics
  }
  
  message("--- Benchmarking Complete ---")
  
  final_results <- dplyr::bind_rows(results_list)
  return(final_results)
}