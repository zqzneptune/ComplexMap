# ComplexMap 2.0.0 (2026-02-19)

## Major Changes (Physical-First Architecture)

*   **[Core Logic]** Shifted to a **Physical-First** landscape philosophy. The network layout now defaults to `alpha = 0.75`, giving 75% weight to physical protein composition to preserve biological specificity and sub-complex architecture.
*   **[Interactive Layer]** Replaced all static plotting functions (`visualizeMapDirectLabels`, `visualizeMapWithLegend`) with a high-performance interactive engine based on **Cytoscape.js** and **Shiny**.
*   **[Visualization]** Added `explore()` and `runComplexMapApp()` to launch the interactive viewer.
*   **[S3 Object]** Overhauled the `print.ComplexMap` method to provide a "Systems Biology Dashboard" highlighting Functional Diversity, Annotation Coverage, and Physical Density.
*   **[Clustering]** Updated functional clustering in `generateNodeAttributes` to use **Average Linkage** and a dynamic scaling factor to force functional diversity and prevent "monochromatic" maps.

## Improvements

*   **[Export]** Enhanced `toCytoscapeJSON` with coordinate scaling and edge pruning to optimize browser performance for large networks.
*   **[Export]** `exportNetwork` now supports GraphML and sanitizes R-specific data types for Cytoscape Desktop compatibility.
*   **[Documentation]** Fully documented missing arguments (`layout_seed`, `...`) and updated vignettes to reflect the v2.0 API.

## Bug Fixes

*   **[Metrics]** Fixed a critical indexing error in `evaluateComplexes` where the Maximum Matching Ratio (MMR) would return `NA` when predicted and reference lists differed in size.
*   **[CRAN]** Eliminated non-ASCII characters in `complexmap_S3.R` and `qcComplexList.R` to ensure cross-platform portability.
*   **[Architecture]** Fixed a corrupt lazy-load database issue by ensuring consistent namespace exports.

---

# ComplexMap 1.1.2 (2026-02-11)

## Added
*   **[createComplexMap]** Added `ifRefineCpx` parameter (default: `FALSE`). Step 1 (Complex Refinement) is now conditional, allowing users to process raw input lists without automatic merging.

## Changed
*   **[generateNodeAttributes]** Enforced a strict 1:1 mapping between Primary Function and Node Color. Removed the weighted color blending logic to ensure visual consistency across the functional landscape.
*   **[Logging]** Updated verbose messaging to explicitly state refinement status during the workflow.

---

# ComplexMap 1.1.1 (2025-12-02)

## Major Improvements
*   **[Parallelism]** Refactored `evaluateComplexes` to use `future.apply` instead of `parallel::mclapply`, enabling reliable multi-core computation on Windows.
*   **[Workflow]** Improved `summarizeThemes()` to accept `add_to_object = TRUE`, allowing theme metadata to be stored directly within the S3 object.

## Bug Fixes
*   **[Dependencies]** Moved `clue` from `Suggests` to `Imports` to ensure the Hungarian algorithm is always available for MMR calculation.
*   **[UI]** Fixed a runtime error in the interactive viewer when the `color.by` argument was null.
*   **[Vignettes]** Resolved build failures in `vignette_02-workflow.Rmd` caused by incorrect internal function calls.

---

# ComplexMap 1.1.0 (2025-11-15)

## Added
*   **[Benchmarking]** Introduced `benchmarkParameters()` to systematically optimize `mergeThreshold` values against reference complexes.
*   **[Metrics]** Added `dice` and `matching_score` similarity options to `refineComplexList`.
*   **[Quantitative]** Added support for a `color.by` argument in visualization functions to map continuous experimental data to node colors.

## Changed
*   **[Defaults]** Updated the default `similarityMethod` in `refineComplexList` to `"matching_score"`.
*   **[Traceability]** `refineComplexList` now returns a `mergeMap` tibble for reproducibility.
*   **[QC]** Overhauled `qcComplexList` output to provide a formatted diagnostic report instead of raw R output.
*   **[Terminology]** Renamed the `"overlap"` similarity metric to `"simpson"` for scientific accuracy.

---

# ComplexMap 1.0.0 (2025-10-10)

*   **Initial stable release.**
*   Introduced the `ComplexMap` S3 object for structured proteomics workflows.
*   Added the high-level `createComplexMap()` wrapper.
*   Implemented core analysis suite: `summarizeThemes()`, `queryMap()`, and `exportNetwork()`.
*   Released comprehensive documentation and three introductory vignettes.