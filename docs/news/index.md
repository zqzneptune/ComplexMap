# Changelog

## ComplexMap 1.1.1 - 2025-12-02

This is a maintenance and stability release that addresses several
critical bugs, improves cross-platform compatibility, and streamlines
the advanced analysis workflow.

### Workflow Enhancements

- The
  [`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
  function has been significantly improved for a more intuitive user
  experience. It now accepts an `add_to_object` argument (defaulting to
  `TRUE`) that adds theme information (`themeId`, `themeLabel`,
  `themePurity`) directly to the `ComplexMap` object’s node table. This
  makes querying by theme using `queryMap(..., type = "theme")` seamless
  and removes the need for manual data merging.

### Major Improvements

- **Cross-Platform Parallel Processing:** The `evaluateComplexes`
  function has been refactored to use the `future.apply` backend,
  replacing the Unix-specific
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html). This
  ensures that parallel computation now works reliably on all operating
  systems, including Windows.

### Bug Fixes

- Fixed a critical runtime error in `visualizeMapInteractive` that
  caused the function to fail when the `color.by` argument was not
  specified. The function now correctly generates tooltips in all use
  cases.

- Resolved a package build failure caused by an incorrect function call
  in the main workflow vignette (`vignette_02-workflow.Rmd`). The
  package can now be built and checked successfully.

- Corrected a dependency issue where `evaluateComplexes` would fail if
  the `clue` package was not installed. `clue` has been moved from
  `Suggests` to `Imports` to guarantee its availability for this core
  function.

### Documentation

- Updated the `README.md` to accurately reflect that `"jaccard"` is the
  default `similarityMethod`, aligning the documentation with the
  package’s diversity-first refinement strategy.

- The advanced analysis vignette (`vignette_03-advanced-analysis.Rmd`)
  has been updated to demonstrate the new, streamlined workflow for
  theme analysis using the refactored
  [`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
  function.

- Synchronized versioning and dates across all documentation files
  (`DESCRIPTION`, `README.md`) for consistency.

## ComplexMap 1.1.0 - 2025-11-15

### Added

- **Quantitative Data Visualization:** All three visualization functions
  (`visualizeMapDirectLabels`, `visualizeMapWithLegend`,
  `visualizeMapInteractive`) now support a `color.by` argument to map
  continuous experimental data directly onto node colors.

- **Robust Parameter Benchmarking:** Added a new
  [`benchmarkParameters()`](https://zqzneptune.github.io/ComplexMap/reference/benchmarkParameters.md)
  function to systematically test and optimize `mergeThreshold` values
  against a reference complex set.

- **New Similarity Metrics:** Added `dice` and `matching_score`
  similarity options to `refineComplexList` for more flexible complex
  merging.

- **New Vignette:** Added a comprehensive vignette
  (`04-quantitative-visualization.Rmd`) to demonstrate the new
  quantitative data workflow.

### Changed

- **New Default Metric:** The default `similarityMethod` in
  `refineComplexList` is now `"matching_score"` to better align
  refinement logic with the standard MMR evaluation metric.

- **Enhanced Traceability:** `refineComplexList` now returns a
  `mergeMap` tibble, providing a clear and reproducible map from
  original to merged complex IDs.

- **Improved QC Report:** Overhauled the output of `qcComplexList` to be
  a user-friendly, formatted report with clear warnings and visual cues,
  rather than raw R output.

- **Metric Renaming:** Renamed the `"overlap"` similarity metric to
  `"simpson"` for clarity and correctness.

### Fixed

- Resolved `R CMD check` errors in vignettes caused by excessive
  parallel core usage in examples (`nCores` is now limited).

- Fixed data type mismatches in vignettes that occurred when joining
  `igraph` membership objects with `dplyr` tibbles.

## ComplexMap 1.0.0

- Initial stable release.
- Complete refactoring of the codebase to adhere to Bioconductor
  standards.
- Introduction of the `ComplexMap` S3 object for a structured workflow.
- Added high-level wrapper
  [`createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
  for end-to-end analysis.
- Added new functions for analysis and exploration:
  [`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md),
  [`queryMap()`](https://zqzneptune.github.io/ComplexMap/reference/queryMap.md),
  and
  [`exportNetwork()`](https://zqzneptune.github.io/ComplexMap/reference/exportNetwork.md).
- Added three comprehensive vignettes.
- Numerous bug fixes and improvements to robustness.
