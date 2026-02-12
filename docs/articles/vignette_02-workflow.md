# 2. Start a Typical Workflow

## Introduction

Protein complexes are the functional machinery of the cell.
High-throughput experimental methods can identify hundreds of putative
protein complexes in a single experiment. The `ComplexMap` package
provides a comprehensive, end-to-end workflow to process, analyze,
annotate, and visualize such a dataset.

This vignette demonstrates the main analysis workflow using the
high-level
[`createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
wrapper function. We will use a dataset of human soluble protein
complexes, originally published in 2012 [A census of human soluble
protein complexes](https://doi.org/10.1016/j.cell.2012.08.011), to
perform a complete analysis from start to finish. We will then showcase
the downstream analysis functions for interpreting and exploring the
results.

## Step 1: Loading Data

The `ComplexMap` package includes the `demoComplexes` dataset, a list of
622 putative human protein complexes. We will also load an example gene
set (GMT) file for functional annotation.

``` r
# Load the example complex list shipped with the package
utils::data("demoComplexes", package = "ComplexMap")

# Get the path to the example GMT file and load it
gmtPath <- ComplexMap::getExampleGmt()
biocartaGmt <- ComplexMap::getGmtFromFile(gmtPath, verbose = FALSE)
```

## Step 2: Running the Main Workflow

The core of the package is the
[`createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
function. This single wrapper function handles all the essential
processing steps:

1.  **Quality Control & Refinement**: Filters complexes by size and
    merges highly redundant ones (defaulting to Jaccard similarity \>
    0.90).
2.  **Functional Enrichment**: Annotates each complex with biological
    functions, prioritizing specificity.
3.  **Network Construction**: Builds a similarity network based on
    shared proteins (physical) and functions.
4.  **Attribute & Topology Generation**: Calculates node colors, sizes,
    and layout coordinates for visualization.

We can pass parameters to the underlying functions directly into this
wrapper.

``` r
# Run the entire workflow with a single command
complexMapObject <- ComplexMap::createComplexMap(
  complexList = demoComplexes,
  gmt = biocartaGmt,
  mergeThreshold = 0.90, # Strict merging to preserve variants
  alpha = 0.75,          # Prioritize physical composition in layout
  verbose = TRUE
)
#> --- Starting ComplexMap Workflow ---
#> Parameters: similarity='jaccard', alpha=0.75 (Diversity Priority)
#> Refinement: Disabled (Using raw input)
#> 
#> Step 1: Skipping refinement (using original complex list)...
#> 
#> Step 2: Running enrichment analysis...
#> Running enrichment for 622 complexes (Cutoff: 0.05)...
#> Annotation complete. Found significant terms for 264 complexes.
#> 
#> Step 3: Building complex network...
#> Building complex network (combined mode, alpha=0.75)...
#> Processing 193131 pairs using 7 cores...
#> Network built: 1681 edges retained.
#> 
#> Step 4: Generating node attributes...
#> Generating node attributes (prioritizing functional specificity)...
#>     -> Clustering 239 terms using co-occurrence (jaccard)
#> Metric: 'jaccard' with unit: 'log'; comparing: 239 vectors
#>     -> Generating diverse palette: 25 functional domains (Average Linkage).
#> 
#> Step 5: Computing map topology...
#> Computing map topology (layout and centrality)...
#> Topology computation complete.
#> 
#> --- ComplexMap Workflow Complete ---
```

The output, `complexMapObject`, is a formal `ComplexMap` S3 object that
contains the complete results of the analysis. Printing the object gives
a high-level systems biology dashboard.

``` r
# Print the object to see a summary
complexMapObject
#> # ComplexMap Object (Physical-First Layout)
#> # ── Physical Structure: 622 nodes, 1681 edges (2.70 edges/node)
#> # ── Functional Landscape:
#> #    • Diversity: 25 distinct functional domains (colors)
#> #    • Coverage:  42.4% of complexes annotated
#> # ── Accessors: `getNodeTable()`, `getEdgeTable()`
#> # ── Analysis:  `summarizeThemes()` to identify physical machines.
```

## Step 3: Summarizing Biological Themes

A key goal of the map is to identify major biological themes. The
[`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
function uses community detection algorithms to find densely connected
network modules (physical neighborhoods) and provides a summary.

``` r
# To get the summary table, we set add_to_object = FALSE
themeSummary <- ComplexMap::summarizeThemes(
  complexMapObject,
  verbose = FALSE,
  add_to_object = FALSE
)

# Display the top 10 largest themes
themeSummary %>%
  dplyr::arrange(dplyr::desc(nodeCount)) %>%
  utils::head(10) %>%
  knitr::kable()
```

| themeId | themeLabel | themePurity | nodeCount | edgeCount |
|---:|:---|---:|---:|---:|
| 2 | BIOCARTA_SALMONELLA_PATHWAY / BIOCARTA_DNAFRAGMENT_PATHWAY | 0.35 | 60 | 167 |
| 1 | BIOCARTA_SALMONELLA_PATHWAY | 0.93 | 38 | 128 |
| 4 | BIOCARTA_MALATEX_PATHWAY | 0.93 | 38 | 58 |
| 6 | BIOCARTA_PELP1_PATHWAY | 0.52 | 36 | 87 |
| 5 | BIOCARTA_SAM68_PATHWAY | 0.55 | 35 | 71 |
| 7 | BIOCARTA_AMAN_PATHWAY | 0.63 | 31 | 101 |
| 9 | BIOCARTA_SALMONELLA_PATHWAY | 0.55 | 29 | 78 |
| 10 | BIOCARTA_RAN_PATHWAY | 0.64 | 29 | 110 |
| 13 | BIOCARTA_PELP1_PATHWAY | 0.60 | 23 | 42 |
| 16 | BIOCARTA_SALMONELLA_PATHWAY | 1.00 | 22 | 24 |

The result is a table listing each theme, its descriptive label (derived
from the most frequent functions within the theme), and its size in
terms of nodes and edges.

## Step 4: Exploring and Querying the Map

The
[`queryMap()`](https://zqzneptune.github.io/ComplexMap/reference/queryMap.md)
function provides a powerful way to programmatically explore the
results.

#### 4.1 Querying for a Specific Protein

Let’s find all complexes that contain the protein “SMAD4”.

``` r
# To ensure our example is robust, let's find a protein to query
# that is guaranteed to be in our final, refined map.
nodes <- ComplexMap::getNodeTable(complexMapObject)
if (nrow(nodes) > 0) {
  first_protein_list <- nodes$proteins
  query_protein <- strsplit(first_protein_list, ",")[[1]][1]
  
  message("Dynamically querying for a protein found in the map: ", query_protein)
  
  protein_complexes <- ComplexMap::queryMap(
    complexMapObject,
    query = query_protein,
    type = "protein"
  )
  
  # Show the primary functional domain of the resulting complexes
  protein_complexes %>%
    dplyr::select(complexId, primaryFunctionalDomain, proteins) %>%
    knitr::kable()
}
#> Dynamically querying for a protein found in the map: PIK3CB
```

| complexId | primaryFunctionalDomain     | proteins      |
|:----------|:----------------------------|:--------------|
| C_17      | BIOCARTA_SALMONELLA_PATHWAY | PIK3CB,PIK3R1 |

#### 4.2 Querying for a Specific Complex

We can also retrieve the data for a single complex of interest.

``` r
if (nrow(nodes) > 0) {
  # Query for the first node in the table
  first_complex_id <- nodes$complexId[1]
  
  complex_data <- ComplexMap::queryMap(
    complexMapObject,
    query = first_complex_id,
    type = "complex"
  )
  
  dplyr::glimpse(complex_data)
}
#> Rows: 1
#> Columns: 12
#> $ complexId               <chr> "C_17"
#> $ proteinCount            <int> 2
#> $ proteins                <chr> "PIK3CB,PIK3R1"
#> $ topIdx                  <int> 6
#> $ primaryFunctionalDomain <chr> "BIOCARTA_SALMONELLA_PATHWAY"
#> $ colorHex                <chr> "#DC716C"
#> $ topEnrichedFunctions    <chr> "BIOCARTA_CDC42RAC_PATHWAY; BIOCARTA_PLC_PATHW…
#> $ sizeMapping             <dbl> 1
#> $ x                       <dbl> -2.053842
#> $ y                       <dbl> 5.333202
#> $ betweenness             <dbl> 0.2034076
#> $ degree                  <dbl> 69
```

## Step 5: Visualization

Finally, we can visualize the entire functional landscape. `ComplexMap`
provides three visualization functions that all work directly with the
`ComplexMap` object.

First, we extract the final node and edge tables for plotting.

``` r
mapLayout <- ComplexMap::getNodeTable(complexMapObject)
networkEdges <- ComplexMap::getEdgeTable(complexMapObject)
```

#### 5.1 Static Plot with a Legend

This version is useful for a clean overview, using a discrete color
legend to represent the functional domains.

``` r
if (nrow(mapLayout) > 0) {
  ComplexMap::visualizeMapWithLegend(mapLayout, networkEdges)
}
#> Visualizing diverse landscape with color legend...
```

![A network plot of protein complexes where node color corresponds to a
functional domain listed in a
legend.](vignette_02-workflow_files/figure-html/vis-legend-1.png)

#### 5.2 Interactive Plot

For deep exploration, an interactive HTML widget is ideal. You can zoom,
pan, and hover over nodes to see detailed tooltips.

``` r
# visNetwork is required for this plot
if (requireNamespace("visNetwork", quietly = TRUE) && nrow(mapLayout) > 0) {
  ComplexMap::visualizeMapInteractive(mapLayout, networkEdges)
}
#> Generating interactive visNetwork plot...
```

## Conclusion

By using the main
[`createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
wrapper, a complete analysis can be run in a single step. The resulting
object can then be easily interpreted, queried, and visualized,
providing a user-friendly and powerful platform for exploring the
landscape of protein complexes.
