# ComplexMap ![ComplexMap hex logo](reference/figures/ComplexMap_logo.png)

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**ComplexMap** is a systems biology R package for the analysis of
protein complex landscapes. It transforms raw lists of putative
complexes into interactive, functionally annotated maps that preserve
the delicate balance between physical composition and biological
hierarchy.

## Core Features (v2.0)

- **Physical-First Architecture:** Layouts are driven primarily by
  **Physical Composition** (shared proteins), ensuring that
  sub-complexes and variants are preserved rather than collapsed into
  generic functional blobs.
- **Interactive Explorer:** Built-in **Cytoscape.js** engine via Shiny.
  Explore your landscape with high-performance zooming, filtering, and
  real-time node inspection.
- **Diversity-First Refinement:** Defaults to conservative **Jaccard**
  similarity (threshold 0.90) to merge only true synonyms while keeping
  distinct biological states separate.
- **Specificity-Driven Annotation:** Automatically groups functional
  terms into diverse domains and calculates **Specificity Scores**
  (`-log10(p) * log2(Fold)`) to prioritize unique labels over generic
  terms.
- **Systems Biology Dashboard:** S3 print methods provide immediate
  feedback on the health of your map (Functional Diversity vs. Physical
  Density).
- **Robust Evaluation:** Compare predictions against gold standards
  (e.g., CORUM) using **Maximum Matching Ratio (MMR)**, PPV, and
  Sensitivity.

------------------------------------------------------------------------

## Installation

You can install the stable version of **ComplexMap** from GitHub.

``` r

# install.packages("remotes")
remotes::install_github("zqzneptune/ComplexMap")
```

------------------------------------------------------------------------

## Workflow at a Glance

This example demonstrates the complete workflow from raw complexes to an
interactive map.

``` r

library(ComplexMap)

# 1. Load demo complexes and example GMT
data("demoComplexes", package = "ComplexMap")
gmt_file <- getExampleGmt()
gmt <- getGmtFromFile(gmt_file, verbose = FALSE)

# 2. Build the ComplexMap
# alpha = 0.75 prioritizes physical structure in the layout
cmap <- createComplexMap(
  complexList = demoComplexes[1:250], 
  gmt = gmt,
  ifRefineCpx = TRUE,
  alpha = 0.75,
  verbose = TRUE
)

# 3. Inspect the Dashboard
cmap
#> # ComplexMap Object (Physical-First Layout)
#> # -- Physical Structure: 224 nodes, 412 edges (1.84 edges/node)
#> # -- Functional Landscape:
#> #    * Diversity: 25 distinct functional domains (colors)
#> #    * Coverage:  42.5% of complexes annotated
#> # -- Accessors: `getNodeTable()`, `getEdgeTable()`
#> # -- Explore:   `explore(cmap)` to launch interactive viewer.

# 4. Explore Interactively
explore(cmap)
```

------------------------------------------------------------------------

## Targeted Exploration

### Powerful Querying

The `queryMap` function supports exact matches for Complex IDs or regex
searches for protein members.

``` r

# Find all complexes containing the protein SMAD4
res <- queryMap(cmap, query = "SMAD4", type = "protein")

# View the dominant functions and membership
dplyr::select(res, complexId, primaryFunctionalDomain, proteinCount)
```

### External Export

Export your landscape to standard formats for high-resolution
publication figures in **Cytoscape Desktop** or **Gephi**.

``` r

# Exports node and edge TSVs sanitized for Cytoscape Desktop
exportNetwork(cmap, filePrefix = "human_proteome_map", format = "cytoscape")

# Export as GraphML
exportNetwork(cmap, filePrefix = "human_proteome_map", format = "graphml")
```

## Performance Benchmarking

Optimize your refinement strategy by testing multiple thresholds against
a reference set (e.g., CORUM).

``` r

data("referenceComplexes", package = "ComplexMap")

bench <- benchmarkParameters(
  complexList = demoComplexes[1:100],
  referenceComplexes = referenceComplexes,
  threshold_range = seq(0.7, 1.0, by = 0.05)
)
```

## Citation

If you use ComplexMap in your research, please cite:

Qingzhou Zhang (2026). **ComplexMap: A Toolset for the Functional
Analysis and Interactive Visualization of Protein Complex Landscapes.**
R package version 2.0.0. <https://github.com/zqzneptune/ComplexMap>

## Contributing

The ComplexMap project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/code_of_conduct/code_of_conduct.md).
We welcome bug reports, feature requests, and pull requests via GitHub
issues.
