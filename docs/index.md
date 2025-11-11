# ComplexMap ![ComplexMap hex logo](reference/figures/ComplexMap_logo.png)

`ComplexMap` provides a comprehensive workflow for the quality control,
refinement, functional enrichment, and network-based analysis of protein
complex datasets. It allows for the quantitative evaluation of predicted
complexes against a reference and provides powerful static and
interactive visualization methods to explore the functional landscape of
a protein complex map.

This package is designed to guide a user from a raw list of putative
protein complexes to a functionally annotated, publication-quality
network visualization.

## Installation

You can install the development version of `ComplexMap` from GitHub.

### 1. Prerequisite Dependencies

Before installing `ComplexMap`, you must first install its dependencies
from CRAN and Bioconductor.

``` r

# Install dependencies from CRAN
install.packages(c(
    "dplyr", "tibble", "magrittr", "rlang", "stringr", "Matrix", 
    "igraph", "ggraph", "tidygraph", "ggplot2", "ggrepel", 
    "RColorBrewer", "scales", "clue", "future", "future.apply", 
    "philentropy", "visNetwork", "devtools"
))

# Install dependencies from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "msigdbr", "AnnotationDbi", "GO.db", "reactome.db", 
    "org.Hs.eg.db" # Example organism database
))
```

### 2. Main Package Installation

Once the dependencies are installed, you can install `ComplexMap`
directly from the `zqzneptune` GitHub repository.

``` r

devtools::install_github("zqzneptune/ComplexMap")
```

## A Quick Example

Here is a minimal example showing the core workflow: from a list of
complexes to a functional network map.

``` r

library(ComplexMap)

# 1. Load the example complex list and a GMT for enrichment
data("demoComplexes")
gmtPath <- getExampleGmt()
gmt <- getGmtFromFile(gmtPath)

# 2. Refine the list (merge redundant complexes)
refinedComplexes <- refineComplexList(demoComplexes, mergeThreshold = 0.75, verbose = FALSE)

# 3. Perform functional enrichment
enrichments <- runComplexEnrichment(refinedComplexes, gmt, verbose = FALSE)

# 4. Build the network
networkEdges <- buildComplexNetwork(refinedComplexes, enrichments, verbose = FALSE)

# 5. Generate node attributes and layout for plotting
nodeAttributes <- generateNodeAttributes(refinedComplexes, enrichments, verbose = FALSE)
mapLayout <- computeMapTopology(nodeAttributes, networkEdges, verbose = FALSE)

# 6. Visualize the final map with a legend

# visualizeMapWithLegend(mapLayout, networkEdges, verbose = FALSE)
```

For a more detailed walkthrough, please see the package vignette:
`browseVignettes("ComplexMap")`.

For a more detailed walkthrough, you can build the vignette while
installing:

``` r

# First, make sure devtools is installed
install.packages("devtools")

# Then install the GitHub package with vignettes
devtools::install_github("zqzneptune/ComplexMap", build_vignettes = TRUE)
```

Then check out the package vignette: `browseVignettes("ComplexMap")`.
