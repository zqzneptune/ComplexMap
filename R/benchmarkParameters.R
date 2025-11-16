#' Benchmark Complex Refinement Parameters
#'
#' @description
#' Evaluates the performance of the complex refinement process across a range of
#' `mergeThreshold` values, providing key metrics against a reference set.
#'
#' @details
#' This function systematically tests different `mergeThreshold` values to help
#' users optimize the complex merging step. For each threshold in `threshold_range`:
#'
#' 1.  It calls `refineComplexList()` with the given threshold.
#' 
#' 2.  It checks the resulting list of refined complexes. If the list is empty
#'     (which can happen with a threshold of 1.0), it records `NA` for all
#'     performance metrics for that threshold.
#'     
#' 3.  If the list is not empty, it calls `evaluateComplexes()` to calculate
#'     PPV, Sensitivity, Accuracy, and MMR.
#'
#' The function returns a single tidy data frame, making it easy to plot the
#' results (e.g., the F1-score) to identify the optimal `mergeThreshold`.
#'
#' @param complexList A list of predicted protein complexes to be refined.
#' @param referenceComplexes A list of reference (gold standard) complexes.
#' @param threshold_range A numeric vector of `mergeThreshold` values to test.
#'   Defaults to `seq(0.5, 1.0, by = 0.05)`.
#' @param ... Additional arguments to be passed down to `refineComplexList`,
#'   such as `similarityMethod` (e.g., "jaccard" or "matching_score"), `minSize`,
#'   or `maxSize`.
#'
#' @return
#' A tidy `tibble` with columns: `mergeThreshold`, `PPV`, `Sn`, `Acc`, and `MMR`.
#' Rows corresponding to thresholds that yielded no complexes will have `NA`
#' values for the metric columns.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # Load the package's demo data
#' data(demoComplexes)
#' data(referenceComplexes)
#'
#' \dontrun{
#' # Run benchmarking over a range of thresholds using the "matching_score" method
#' benchmark_results <- benchmarkParameters(
#'   complexList = demoComplexes,
#'   referenceComplexes = referenceComplexes,
#'   threshold_range = seq(0.7, 1.0, by = 0.05),
#'   similarityMethod = "matching_score"
#' )
#'
#' print(benchmark_results)
#'
#' # To find the optimal threshold, calculate the F1-score (harmonic mean
#' # of PPV and Sn) and plot the results.
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE)) {
#'
#'   plot_data <- benchmark_results %>%
#'     dplyr::mutate(F1_Score = 2 * (PPV * Sn) / (PPV + Sn))
#'
#'   # Find the best threshold
#'   best_threshold <- plot_data %>%
#'     dplyr::filter(F1_Score == max(F1_Score, na.rm = TRUE)) %>%
#'     dplyr::pull(mergeThreshold)
#'
#'   ggplot2::ggplot(plot_data, ggplot2::aes(x = mergeThreshold, y = F1_Score)) +
#'     ggplot2::geom_line(color = "blue") +
#'     ggplot2::geom_point(size = 3, color = "blue") +
#'     ggplot2::geom_vline(
#'       xintercept = best_threshold, linetype = "dashed", color = "red"
#'     ) +
#'     ggplot2::scale_x_continuous(breaks = plot_data$mergeThreshold) +
#'     ggplot2::labs(
#'       title = "F1-Score vs. Merge Threshold",
#'       subtitle = paste("Optimal threshold based on F1-Score:", best_threshold),
#'       x = "Merge Threshold",
#'       y = "F1-Score"
#'     ) +
#'     ggplot2::theme_minimal(base_size = 14)
#' }
#' }
#'
benchmarkParameters <- function(complexList, referenceComplexes,
                                threshold_range = seq(0.5, 1.0, by = 0.05),
                                ...) {

  # Initialize a list to store the results for each iteration
  results_list <- vector("list", length(threshold_range))

  message("--- Starting Parameter Benchmarking ---")

  # Loop over each specified threshold value
  for (i in seq_along(threshold_range)) {
    current_threshold <- threshold_range[i]
    if (interactive()) {
        message(sprintf("Benchmarking with mergeThreshold = %.2f...",
                        current_threshold))
    }


    # Call refineComplexList, passing through any additional arguments
    # from '...' (like similarityMethod). Verbose is turned off for the loop.
    refinement_output <- refineComplexList(
      complexList = complexList,
      mergeThreshold = current_threshold,
      verbose = FALSE,
      ...
    )

    refined_complexes <- refinement_output$refinedComplexes

    # Check if the refinement process yielded any complexes
    if (length(refined_complexes) == 0) {
      # If the list is empty, record NA for all metrics
      metrics <- tibble::tibble(
        mergeThreshold = current_threshold,
        PPV = NA_real_,
        Sn = NA_real_,
        Acc = NA_real_,
        MMR = NA_real_
      )
    } else {
      # If complexes remain, evaluate them against the reference
      eval_metrics <- evaluateComplexes(
        predictedComplexes = refined_complexes,
        referenceComplexes = referenceComplexes,
        verbose = FALSE # Suppress verbose output within the loop
      )

      # Combine the threshold and the metrics into a single tibble row
      metrics <- tibble::tibble(
        mergeThreshold = current_threshold,
        PPV = eval_metrics$PPV,
        Sn = eval_metrics$Sn,
        Acc = eval_metrics$Acc,
        MMR = eval_metrics$MMR
      )
    }
    results_list[[i]] <- metrics
  }

  message("--- Benchmarking Complete ---")

  # Combine all the individual results into a single final tibble
  final_results <- dplyr::bind_rows(results_list)
  return(final_results)
}