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
    merges highly redundant ones.

2.  **Functional Enrichment**: Annotates each complex with biological
    functions.

3.  **Network Construction**: Builds a similarity network based on
    shared proteins and functions.

4.  **Attribute & Topology Generation**: Calculates node colors, sizes,
    and layout coordinates for visualization.

We can pass parameters to the underlying functions (like
`mergeThreshold` for `refineComplexList`) directly into this wrapper.

``` r
# Run the entire workflow with a single command
# We will merge complexes with a Jaccard similarity of 0.75 or higher
complexMapObject <- ComplexMap::createComplexMap(
  complexList = demoComplexes,
  gmt = biocartaGmt,
  mergeThreshold = 0.75
)
#> --- Starting ComplexMap Workflow ---
#> 
#> Step 1: Refining complex list...
#> 
#> --- Refining Input Complex List ---
#> Filtered 112 complexes by size. Retaining 510.
#> Identifying merge groups with matching_score >= 0.75...
#> Found 0 merge groups. Merging 0 complexes into 510.
#> Merging complete. Final list has 510 complexes.
#> 
#> --- Refinement Complete ---
#> 
#> Step 2: Running enrichment analysis...
#> Running enrichment for 510 complexes...
#> Annotation complete. Found terms for 239 complexes.
#> 
#> Step 3: Building complex network...
#> Building complex network using 'jaccard' similarity...
#> Using 7 cores for parallel processing.
#> Processing 129795 complex pairs...
#> Split into 130 chunks of up to 1000 pairs each.
#> Combining results from chunks...
#> Calculating final weights and filtering...
#> Network construction complete: 1390 edges retained.
#> 
#> Step 4: Generating node attributes...
#> Generating core node attributes (function and color)...
#>     -> Clustering 231 terms using co-occurrence (overlap)
#>     -> Identified 3 functional domains
#> 
#> Step 5: Computing map topology...
#> Computing map topology (layout and centrality)...
#> Topology computation complete.
#> 
#> --- ComplexMap Workflow Complete ---
```

The output, `complexMapObject`, is a formal `ComplexMap` S3 object that
contains the complete results of the analysis. Printing the object gives
a high-level summary.

``` r
# Print the object to see a summary
complexMapObject
#> # A ComplexMap Object
#> # ── 510 nodes and 1390 edges
#> # ── 3 major biological themes identified.
#> # ── Use `getNodeTable()` or `getEdgeTable()` to access data.
```

## Step 3: Summarizing Biological Themes

A key goal of the map is to identify major biological themes. The
[`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
function uses community detection algorithms to find densely connected
network modules and provides a summary.

``` r
themeSummary <- ComplexMap::summarizeThemes(complexMapObject)
#> Summarizing themes using the 'louvain' community algorithm...
#> Identified 98 distinct themes.

# Display the top 10 largest themes
themeSummary %>%
  dplyr::arrange(dplyr::desc(nodeCount)) %>%
  utils::head(10)
#> # A tibble: 10 × 4
#>    themeId themeLabel                  nodeCount edgeCount
#>      <int> <chr>                           <int>     <int>
#>  1       1 BIOCARTA_PTDINS_PATHWAY            56       180
#>  2       2 BIOCARTA_PROTEASOME_PATHWAY        47        86
#>  3       8 BIOCARTA_PROTEASOME_PATHWAY        45        68
#>  4       9 BIOCARTA_PROTEASOME_PATHWAY        31        68
#>  5       6 BIOCARTA_PROTEASOME_PATHWAY        27        47
#>  6       5 BIOCARTA_PROTEASOME_PATHWAY        25        78
#>  7      11 BIOCARTA_PROTEASOME_PATHWAY        25        92
#>  8      12 BIOCARTA_PROTEASOME_PATHWAY        24       101
#>  9       4 BIOCARTA_PROTEASOME_PATHWAY        20        49
#> 10      10 BIOCARTA_CACAM_PATHWAY             20        20
```

The result is a table listing each theme, its descriptive label (derived
from the most common function within the theme), and its size in terms
of nodes and edges.

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
first_protein_list <- nodes$proteins[1]
query_protein <- strsplit(first_protein_list, ",")[[1]][1]

message("Dynamically querying for a protein found in the map: ", query_protein)
#> Dynamically querying for a protein found in the map: PRKACB

protein_complexes <- ComplexMap::queryMap(
  complexMapObject,
  query = query_protein,
  type = "protein"
)

# Show the primary functional domain of the resulting complexes
protein_complexes %>%
  dplyr::select(complexId, primaryFunctionalDomain, proteins)
#> # A tibble: 1 × 3
#>   complexId   primaryFunctionalDomain proteins                           
#>   <chr>       <chr>                   <chr>                              
#> 1 CpxMap_0359 BIOCARTA_CACAM_PATHWAY  PRKACB,PRKACA,PRKAR2A,CIRBP,CAPRIN1
```

#### 4.2 Querying for a Specific Complex

We can also retrieve the data for a single complex of interest.

``` r
# Note: The exact CpxMap ID may vary slightly between runs
# if refinement results change. We query for the first node in the table.
first_complex_id <- ComplexMap::getNodeTable(complexMapObject)$complexId[1]

complex_data <- ComplexMap::queryMap(
  complexMapObject,
  query = first_complex_id,
  type = "complex"
)

dplyr::glimpse(complex_data)
#> Rows: 1
#> Columns: 11
#> $ complexId               <chr> "CpxMap_0359"
#> $ proteinCount            <int> 5
#> $ proteins                <chr> "PRKACB,PRKACA,PRKAR2A,CIRBP,CAPRIN1"
#> $ primaryFunctionalDomain <chr> "BIOCARTA_CACAM_PATHWAY"
#> $ topEnrichedFunctions    <chr> "BIOCARTA_AGPCR_PATHWAY; BIOCARTA_AKAP13_PATHW…
#> $ colorHex                <chr> "#836FA8"
#> $ sizeMapping             <dbl> 2.321928
#> $ x                       <dbl> 7.989894
#> $ y                       <dbl> 6.128077
#> $ betweenness             <dbl> 0.1544715
#> $ degree                  <dbl> 33
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
ComplexMap::visualizeMapWithLegend(mapLayout, networkEdges)
#> Visualizing ComplexMap with a color legend...
```

![A network plot of protein complexes where node color corresponds to a
functional domain listed in a
legend.](vignette_02-workflow_files/figure-html/vis-legend-1.png)

#### 5.2 Interactive Plot

For deep exploration, an interactive HTML widget is ideal. You can zoom,
pan, and hover over nodes to see detailed tooltips.

``` r
# visNetwork is required for this plot
if (requireNamespace("visNetwork", quietly = TRUE)) {
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
