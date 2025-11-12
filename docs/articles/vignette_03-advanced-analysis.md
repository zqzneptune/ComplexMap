# 3. Explore Advanced Analyses

## Introduction

While the
[`ComplexMap::createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
function provides a seamless, one-step workflow for most use cases, the
package also includes a suite of granular functions for users who
require more control, wish to benchmark their results, or need to export
data for other tools.

This vignette covers these advanced topics:

1.  **Granular Workflow Control**: Running the quality control and
    refinement steps manually to inspect intermediate results.

2.  **Benchmarking Predictions**: Using
    [`evaluateComplexes()`](https://zqzneptune.github.io/ComplexMap/reference/evaluateComplexes.md)
    to compare your experimental complexes against a reference set like
    CORUM.

3.  **Advanced Network Exploration**: Adding theme information back to a
    `ComplexMap` object to enable powerful, theme-based queries.

4.  **Interoperability**: Exporting the final network for visualization
    and analysis in external software like Cytoscape.

We will use the package’s built-in datasets for this demonstration.

``` r

# Load the experimental and reference complex lists
utils::data("demoComplexes", package = "ComplexMap")
utils::data("referenceComplexes", package = "ComplexMap")
```

## 1. Granular Workflow Control

The
[`createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
function automatically runs quality control and refinement. However, you
can run these steps manually using the underlying functions to better
understand and tune your analysis.

#### 1.1 Quality Control

First, we can run
[`qcComplexList()`](https://zqzneptune.github.io/ComplexMap/reference/qcComplexList.md)
on our raw input data to get a detailed report on complex sizes and
redundancy.

``` r

# Run the QC function with verbose output
ComplexMap::qcComplexList(demoComplexes)
#> 
#> --- Running Quality Control on Complex List ---
#> 
#> [1] Basic Statistics:
#>     - Total number of complexes: 622
#>     - Total number of unique proteins: 2649
#> 
#> [2] Complex Size Distribution:
#>        Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#>       1.000   3.000   4.000   5.183   5.000 102.000
#> Warning in ComplexMap::qcComplexList(demoComplexes): 112 complexes have fewer
#> than 3 members.
#> 
#> [3] Redundancy Analysis:
#>     - Distribution of Jaccard similarity scores:
#>          Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#>     0.0000000 0.0000000 0.0000000 0.0003226 0.0000000 0.5000000
#>     - No highly redundant complex pairs detected.
#> 
#> --- QC Complete ---
```

The report provides summary statistics and warns us that a number of
complex pairs are highly redundant (Jaccard similarity \>= 0.8). This
confirms that the refinement step is necessary.

#### 1.2 Manual Refinement

Next, we can manually run
[`refineComplexList()`](https://zqzneptune.github.io/ComplexMap/reference/refineComplexList.md).
This allows us to inspect the set of complexes *after* merging but
*before* functional enrichment, which can be useful for debugging or
parameter tuning.

``` r

message("Number of complexes before refinement: ", length(demoComplexes))
#> Number of complexes before refinement: 622

# Merge complexes with a Jaccard similarity of 0.75 or higher
refinedComplexes <- ComplexMap::refineComplexList(
  demoComplexes,
  mergeThreshold = 0.75,
  verbose = FALSE # Set to TRUE for detailed messages
)

message("Number of complexes after refinement: ", length(refinedComplexes))
#> Number of complexes after refinement: 510
```

Running this step manually gives us direct access to the
`refinedComplexes` list for further inspection before proceeding with
the rest of the workflow.

## 2. Benchmarking Predictions

For methods development, it is essential to evaluate how well a set of
predicted complexes matches a “gold standard” reference. The
[`evaluateComplexes()`](https://zqzneptune.github.io/ComplexMap/reference/evaluateComplexes.md)
function calculates four standard metrics for this purpose.

We will compare our `refinedComplexes` against the `referenceComplexes`
(a curated subset of CORUM).

``` r

# Evaluate the refined complexes against the reference set
evaluationMetrics <- ComplexMap::evaluateComplexes(
  predictedComplexes = refinedComplexes,
  referenceComplexes = referenceComplexes,
  nCores = 2,
  verbose = FALSE
)

# Print the resulting metrics
print(evaluationMetrics)
#> $PPV
#> [1] 0.8005051
#> 
#> $Sn
#> [1] 0.2510583
#> 
#> $Acc
#> [1] 0.4483006
#> 
#> $MMR
#> [1] 0.09403458
```

- **PPV** (Positive Predictive Value): Measures the accuracy of your
  predictions.

- **Sn** (Sensitivity): Measures the coverage of the reference set.

- **Acc** (Accuracy): The geometric mean of PPV and Sn.

- **MMR** (Maximum Matching Ratio): A score based on finding the best
  one-to-one mapping between predicted and reference complexes.

These metrics provide a quantitative assessment of your complex
prediction quality.

## 3. Advanced Network Exploration: Querying by Theme

In the main workflow vignette, we used
[`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
to get a table of biological themes. A powerful advanced technique is to
add this theme information back into the `ComplexMap` object itself,
enabling direct queries.

First, let’s generate a `ComplexMap` object.

``` r

gmtPath <- ComplexMap::getExampleGmt()
gmt <- ComplexMap::getGmtFromFile(gmtPath, verbose = FALSE)

# We use the manually refined list as input
cm_obj <- ComplexMap::createComplexMap(
  complexList = refinedComplexes,
  gmt = gmt,
  verbose = FALSE
)
#> Metric: 'jaccard' with unit: 'log'; comparing: 231 vectors
```

Now, we generate the themes and add the theme assignments to our node
table.

``` r

# Generate the theme summary to get theme labels
theme_summary <- ComplexMap::summarizeThemes(cm_obj, verbose = FALSE)

# To get the theme assignment for each node, we must re-run community detection
nodes <- ComplexMap::getNodeTable(cm_obj)
edges <- ComplexMap::getEdgeTable(cm_obj)
graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
communities <- igraph::cluster_louvain(graph)

# Create a mapping from themeId (integer) to themeLabel (character)
theme_map <- theme_summary[, c("themeId", "themeLabel")]

# Add theme ID and label columns to our node table
nodes$themeId <- igraph::membership(communities)
nodes <- dplyr::left_join(nodes, theme_map, by = "themeId")

# **Crucially, update the object with the new node information**
cm_obj$nodes <- nodes
```

With the prepared object, we can now find all complexes belonging to a
specific theme, for example, the “26S Proteasome”.

``` r

# To make the example robust, we query for a theme label that we know
# exists because we just generated it in the `theme_summary` table.
query_label <- theme_summary$themeLabel

message("Dynamically querying for a theme found in the map: ", query_label)
#> Dynamically querying for a theme found in the map: BIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_SAM68_PATHWAYBIOCARTA_ETS_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EPONFKB_PATHWAYBIOCARTA_MRP_PATHWAYBIOCARTA_EICOSANOID_PATHWAYBIOCARTA_EICOSANOID_PATHWAYMixed or UnenrichedMixed or UnenrichedBIOCARTA_MCM_PATHWAYBIOCARTA_RAB_PATHWAYBIOCARTA_EICOSANOID_PATHWAYMixed or UnenrichedBIOCARTA_NUCLEARRS_PATHWAYBIOCARTA_S1P_PATHWAYMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedBIOCARTA_GABA_PATHWAYMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedBIOCARTA_LDL_PATHWAYMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedMixed or UnenrichedBIOCARTA_ARAP_PATHWAYMixed or UnenrichedMixed or UnenrichedMixed or Unenriched

theme_complexes <- ComplexMap::queryMap(
  cm_obj,
  query = query_label,
  type = "theme"
)
#> Warning in nodes$themeLabel == query: longer object length is not a multiple of
#> shorter object length

# Show the results
theme_complexes %>%
  dplyr::select(complexId, themeLabel, proteinCount, betweenness)
#> # A tibble: 114 × 4
#>    complexId   themeLabel                  proteinCount betweenness
#>    <chr>       <chr>                              <int>       <dbl>
#>  1 CpxMap_0359 BIOCARTA_EICOSANOID_PATHWAY            5      0.154 
#>  2 CpxMap_0414 BIOCARTA_EICOSANOID_PATHWAY            8      0.136 
#>  3 CpxMap_0508 BIOCARTA_EICOSANOID_PATHWAY           42      0.111 
#>  4 CpxMap_0401 BIOCARTA_EICOSANOID_PATHWAY            6      0.100 
#>  5 CpxMap_0509 BIOCARTA_EICOSANOID_PATHWAY          102      0.0918
#>  6 CpxMap_0090 BIOCARTA_EICOSANOID_PATHWAY            4      0.0888
#>  7 CpxMap_0143 BIOCARTA_EICOSANOID_PATHWAY            4      0.0700
#>  8 CpxMap_0182 BIOCARTA_EICOSANOID_PATHWAY            4      0.0657
#>  9 CpxMap_0501 BIOCARTA_EICOSANOID_PATHWAY           20      0.0598
#> 10 CpxMap_0372 BIOCARTA_EICOSANOID_PATHWAY            6      0.0547
#> # ℹ 104 more rows
```

This technique provides a powerful way to programmatically subset and
analyze the biological modules within your map.

## 4. Interoperability: Exporting to Cytoscape

While `ComplexMap` provides excellent internal visualizations, you may
want to use other powerful tools like
[Cytoscape](https://cytoscape.org/) for advanced analysis, network
manipulation, or creating publication-quality figures.

The
[`exportNetwork()`](https://zqzneptune.github.io/ComplexMap/reference/exportNetwork.md)
function writes the node and edge tables to disk as tab-separated files
(`.tsv`), which can be easily imported into Cytoscape.

``` r

# Use a temporary directory for this example
temp_dir <- tempdir()
file_prefix <- file.path(temp_dir, "human_complex_map")

# Export the network
ComplexMap::exportNetwork(cm_obj, filePrefix = file_prefix)
#> Exporting network in Cytoscape format to:
#>  -> /var/folders/x8/ngpdbcfd3pj5sf3r9kkxgs740000gn/T//RtmpZDxIxh/human_complex_map_nodes.tsv
#>  -> /var/folders/x8/ngpdbcfd3pj5sf3r9kkxgs740000gn/T//RtmpZDxIxh/human_complex_map_edges.tsv
#> Export complete.

# List the files that were created
list.files(temp_dir, pattern = "human_complex_map")
#> [1] "human_complex_map_edges.tsv" "human_complex_map_nodes.tsv"
```

To use these in Cytoscape, you would:

1.  Go to `File > Import > Network from File...` and select the
    `_edges.tsv` file.

2.  In the mapping dialog, map the `source_complex_id` column to “Source
    Node” and `target_complex_id` to “Target Node”.

3.  Once the network is loaded, go to
    `File > Import > Table from File...` and select the `_nodes.tsv`
    file.

4.  Ensure the “Key column for network” is set to `complexId`.

This will import all your node attributes (colors, sizes, functional
domains, etc.) and map them directly onto your network.

## Conclusion

This vignette has demonstrated the advanced capabilities of the
`ComplexMap` package. Beyond the simple one-step workflow, the package
provides granular functions for quality control, rigorous benchmarking
tools for methods development, and powerful utilities for deep,
theme-based querying and exporting results to other platforms. This
layered functionality makes `ComplexMap` a comprehensive toolkit for
both novice and expert users.
