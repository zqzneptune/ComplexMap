# Benchmark Complex Refinement Parameters

Evaluates the performance of the complex refinement process across a
range of \`mergeThreshold\` values, providing key metrics against a
reference set.

## Usage

``` r
benchmarkParameters(
  complexList,
  referenceComplexes,
  threshold_range = seq(0.5, 1, by = 0.05),
  ...
)
```

## Arguments

- complexList:

  A list of predicted protein complexes to be refined.

- referenceComplexes:

  A list of reference (gold standard) complexes.

- threshold_range:

  A numeric vector of \`mergeThreshold\` values to test. Defaults to
  \`seq(0.5, 1.0, by = 0.05)\`.

- ...:

  Additional arguments to be passed down to \`refineComplexList\`, such
  as \`similarityMethod\` (e.g., "jaccard" or "overlap"), \`minSize\`,
  or \`maxSize\`.

## Value

A tidy \`tibble\` with columns: \`mergeThreshold\`, \`PPV\`, \`Sn\`,
\`Acc\`, and \`MMR\`. Rows corresponding to thresholds that yielded no
complexes will have \`NA\` values for the metric columns.

## Details

This function systematically tests different \`mergeThreshold\` values
to help users optimize the complex merging step. For each threshold in
\`threshold_range\`:

1\. It calls \`refineComplexList()\` with the given threshold.

2\. It checks the resulting list of refined complexes. If the list is
empty (which can happen with a threshold of 1.0), it records \`NA\` for
all performance metrics for that threshold.

3\. If the list is not empty, it calls \`evaluateComplexes()\` to
calculate PPV, Sensitivity, Accuracy, and MMR.

The function returns a single tidy data frame, making it easy to plot
the results (e.g., the F1-score) to identify the optimal
\`mergeThreshold\`.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Load the package's demo data
data(demoComplexes)
data(referenceComplexes)

if (FALSE) { # \dontrun{
# Run benchmarking over a range of thresholds using the "overlap" method
benchmark_results <- benchmarkParameters(
  complexList = demoComplexes,
  referenceComplexes = referenceComplexes,
  threshold_range = seq(0.7, 1.0, by = 0.05),
  similarityMethod = "overlap"
)

print(benchmark_results)

# To find the optimal threshold, calculate the F1-score (harmonic mean
# of PPV and Sn) and plot the results.
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("dplyr", quietly = TRUE)) {

  plot_data <- benchmark_results %>%
    dplyr::mutate(F1_Score = 2 * (PPV * Sn) / (PPV + Sn))

  # Find the best threshold
  best_threshold <- plot_data %>%
    dplyr::filter(F1_Score == max(F1_Score, na.rm = TRUE)) %>%
    dplyr::pull(mergeThreshold)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = mergeThreshold, y = F1_Score)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_point(size = 3, color = "blue") +
    ggplot2::geom_vline(
      xintercept = best_threshold, linetype = "dashed", color = "red"
    ) +
    ggplot2::scale_x_continuous(breaks = plot_data$mergeThreshold) +
    ggplot2::labs(
      title = "F1-Score vs. Merge Threshold",
      subtitle = paste("Optimal threshold based on F1-Score:", best_threshold),
      x = "Merge Threshold",
      y = "F1-Score"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}
} # }
```
