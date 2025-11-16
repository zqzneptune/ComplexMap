# ComplexMap 1.1.0 - 2025-11-15

## Added

-   **Quantitative Data Visualization:** All three visualization functions (`visualizeMapDirectLabels`, `visualizeMapWithLegend`, `visualizeMapInteractive`) now support a `color.by` argument to map continuous experimental data directly onto node colors.

-   **Robust Parameter Benchmarking:** Added a new `benchmarkParameters()` function to systematically test and optimize `mergeThreshold` values against a reference complex set.

-   **New Similarity Metrics:** Added `dice` and `matching_score` similarity options to `refineComplexList` for more flexible complex merging.

-   **New Vignette:** Added a comprehensive vignette (`04-quantitative-visualization.Rmd`) to demonstrate the new quantitative data workflow.

## Changed

-   **New Default Metric:** The default `similarityMethod` in `refineComplexList` is now `"matching_score"` to better align refinement logic with the standard MMR evaluation metric.

-   **Enhanced Traceability:** `refineComplexList` now returns a `mergeMap` tibble, providing a clear and reproducible map from original to merged complex IDs.

-   **Improved QC Report:** Overhauled the output of `qcComplexList` to be a user-friendly, formatted report with clear warnings and visual cues, rather than raw R output.

-   **Metric Renaming:** Renamed the `"overlap"` similarity metric to `"simpson"` for clarity and correctness.

## Fixed

-   Resolved `R CMD check` errors in vignettes caused by excessive parallel core usage in examples (`nCores` is now limited).

-   Fixed data type mismatches in vignettes that occurred when joining `igraph` membership objects with `dplyr` tibbles.

# ComplexMap 1.0.0

* Initial stable release.
* Complete refactoring of the codebase to adhere to Bioconductor standards.
* Introduction of the `ComplexMap` S3 object for a structured workflow.
* Added high-level wrapper `createComplexMap()` for end-to-end analysis.
* Added new functions for analysis and exploration: `summarizeThemes()`, `queryMap()`, and `exportNetwork()`.
* Added three comprehensive vignettes.
* Numerous bug fixes and improvements to robustness.