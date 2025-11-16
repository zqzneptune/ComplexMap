# 3. Advanced Analysis and Visualization

## Advanced Analysis Workflow

Once you have a `ComplexMap` object, the package provides several
functions for deeper analysis and interpretation of the network
structure. This vignette covers:

1.  **Optimizing Parameters**: Using
    [`evaluateComplexes()`](https://zqzneptune.github.io/ComplexMap/reference/evaluateComplexes.md)
    to tune the `refineComplexList` step.

2.  **Summarizing Themes**: Using
    [`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
    to identify major biological clusters.

3.  **Querying the Map**: Using
    [`queryMap()`](https://zqzneptune.github.io/ComplexMap/reference/queryMap.md)
    to find specific proteins, complexes, or themes.

``` r
library(ComplexMap)
library(dplyr)
library(ggplot2)
library(igraph)
```

### 1. Optimizing Refinement with `evaluateComplexes`

The `mergeThreshold` in `refineComplexList` is a critical parameter that
determines how aggressively smaller complexes are merged. A low
threshold creates fewer, larger complexes, while a high threshold
preserves more distinct complexes. The optimal value often depends on
the dataset.

The
[`evaluateComplexes()`](https://zqzneptune.github.io/ComplexMap/reference/evaluateComplexes.md)
function can be used to score a refined complex list against a gold
standard (like the `referenceComplexes` dataset) to help guide this
choice. It returns four key metrics: Positive Predictive Value (PPV),
Sensitivity (Sn), Accuracy (Acc), and Maximum Matching Ratio (MMR).

Let’s test a `mergeThreshold` of 0.85.

``` r
# Load the package's demo and reference data
data(demoComplexes)
data(referenceComplexes)

# refineComplexList now returns a list with two elements
refinement_output <- refineComplexList(
  demoComplexes,
  mergeThreshold = 0.85,
  verbose = FALSE
)

# We must now explicitly use the `$refinedComplexes` element for evaluation
metrics <- evaluateComplexes(
  predictedComplexes = refinement_output$refinedComplexes,
  referenceComplexes = referenceComplexes,
  nCores = 2, # <-- THIS IS THE FIX: Limit cores for check environment
  verbose = FALSE
)

# Display the resulting metrics
knitr::kable(as.data.frame(metrics))
```

|       PPV |        Sn |       Acc |       MMR |
|----------:|----------:|----------:|----------:|
| 0.8005051 | 0.2510583 | 0.4483006 | 0.0940346 |

For more systematic tuning, see the documentation for the new
[`benchmarkParameters()`](https://zqzneptune.github.io/ComplexMap/reference/benchmarkParameters.md)
function, which automates this process over a range of thresholds.

### 2. Summarizing Biological Themes

A key goal of `ComplexMap` is to organize complexes into a functional
landscape. The
[`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
function uses network community detection algorithms (like Louvain, the
default) to find densely connected clusters of complexes within the map.
It then automatically labels each cluster based on the most common
functional annotation among its members.

Let’s generate a map with a slightly lower merge threshold to create a
more interconnected network for this example.

``` r
# We'll use the example GMT file that comes with the package
gmtPath <- getExampleGmt()
gmt <- getGmtFromFile(gmtPath, verbose = FALSE)

# Create the map object
cm_obj_adv <- createComplexMap(
  complexList = demoComplexes,
  gmt = gmt,
  mergeThreshold = 0.75, # Lower threshold for more connections
  verbose = FALSE
)
```

Now, we can run
[`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
on this object.

``` r
# The function returns a tibble summarizing each detected theme
theme_summary <- summarizeThemes(cm_obj_adv, verbose = FALSE)

# Display the largest themes
theme_summary %>%
  arrange(desc(nodeCount)) %>%
  head(10) %>%
  knitr::kable()
```

| themeId | themeLabel                  | nodeCount | edgeCount |
|--------:|:----------------------------|----------:|----------:|
|       1 | BIOCARTA_PTDINS_PATHWAY     |        56 |       180 |
|       2 | BIOCARTA_PROTEASOME_PATHWAY |        47 |        86 |
|       8 | BIOCARTA_PROTEASOME_PATHWAY |        45 |        68 |
|       9 | BIOCARTA_PROTEASOME_PATHWAY |        31 |        68 |
|       6 | BIOCARTA_PROTEASOME_PATHWAY |        27 |        47 |
|       5 | BIOCARTA_PROTEASOME_PATHWAY |        25 |        78 |
|      11 | BIOCARTA_PROTEASOME_PATHWAY |        25 |        92 |
|      12 | BIOCARTA_PROTEASOME_PATHWAY |        24 |       101 |
|       4 | BIOCARTA_PROTEASOME_PATHWAY |        20 |        49 |
|      10 | BIOCARTA_CACAM_PATHWAY      |        20 |        20 |

This summary shows us the major biological stories in our dataset, such
as the large “Ribosome” and “Proteasome” communities.

### 3. Querying the Map for Specific Information

The
[`queryMap()`](https://zqzneptune.github.io/ComplexMap/reference/queryMap.md)
function provides an easy way to retrieve specific nodes from the
`ComplexMap` object based on different criteria.

#### Query by Protein

Find all complexes that contain a specific protein.

``` r
# Find all complexes containing the protein "UBA1"
queryMap(cm_obj_adv, query = "UBA1", type = "protein") %>%
  select(complexId, primaryFunctionalDomain, proteinCount, proteins) %>%
  knitr::kable()
```

| complexId | primaryFunctionalDomain | proteinCount | proteins |
|:---|:---|---:|:---|
| CpxMap_0497 | BIOCARTA_PROTEASOME_PATHWAY | 18 | UBA1,TRMT1,RARS1,QARS1,NARS1,MARS1,LARS1,KARS1,IARS1,GARS1,EPRS1,DARS1,AARS1,VBP1,EEF1E1,ASNS,AIMP2,AIMP1 |
| CpxMap_0195 | BIOCARTA_PROTEASOME_PATHWAY | 4 | WDR73,USH1C,UBA1,NMD3 |

#### Query by Complex ID

Retrieve the data for a single complex.

``` r
# Get the data for complex "CpxMap_0001"
queryMap(cm_obj_adv, query = "CpxMap_0001", type = "complex") %>%
  knitr::kable()
```

| complexId | proteinCount | proteins | primaryFunctionalDomain | topEnrichedFunctions | colorHex | sizeMapping | x | y | betweenness | degree |
|:---|---:|:---|:---|:---|:---|---:|---:|---:|---:|---:|
| CpxMap_0001 | 3 | WDR82,ULK3,UCHL1 | Unenriched | NA | \#CCCCCC | 1.584963 | 6.421011 | -14.1341 | 0 | 0 |

#### Query by Theme

To query by theme, we first need to add the theme labels from
[`summarizeThemes()`](https://zqzneptune.github.io/ComplexMap/reference/summarizeThemes.md)
back into our main `ComplexMap` object’s node table.

``` r
# To add the theme labels to our main node table for querying or visualization,
# we first need to get the theme ID for each node.

# Re-run the community detection to get the membership vector
graph <- igraph::graph_from_data_frame(
  d = getEdgeTable(cm_obj_adv), 
  vertices = getNodeTable(cm_obj_adv), 
  directed = FALSE
)
communities <- igraph::cluster_louvain(graph)

# Create a tibble mapping complexId to themeId
# Coerce the 'membership' object to a simple integer for joining
node_themes <- tibble(
  complexId = names(igraph::membership(communities)),
  themeId = as.integer(igraph::membership(communities))
)
theme_summary$themeId <-
  as.integer(theme_summary$themeId)

# Join this mapping with the theme summary table to get the labels
node_theme_labels <- node_themes %>%
  left_join(theme_summary, by = "themeId")

# Finally, join the theme labels back to the main node table of our object
nodes_with_themes <- getNodeTable(cm_obj_adv) %>%
  left_join(node_theme_labels, by = "complexId")

# Update the ComplexMap object with this new node table
cm_obj_adv$nodes <- nodes_with_themes

# Preview the new columns
head(select(cm_obj_adv$nodes, complexId, themeId, themeLabel)) %>%
  knitr::kable()
```

| complexId   | themeId | themeLabel                  |
|:------------|--------:|:----------------------------|
| CpxMap_0359 |       1 | BIOCARTA_PTDINS_PATHWAY     |
| CpxMap_0414 |       2 | BIOCARTA_PROTEASOME_PATHWAY |
| CpxMap_0508 |       3 | BIOCARTA_PROTEASOME_PATHWAY |
| CpxMap_0401 |       4 | BIOCARTA_PROTEASOME_PATHWAY |
| CpxMap_0509 |       5 | BIOCARTA_PROTEASOME_PATHWAY |
| CpxMap_0090 |       1 | BIOCARTA_PTDINS_PATHWAY     |

Now that the `themeLabel` column exists in our node table, we can query
it.

``` r
# Find all complexes belonging to the "Proteasome" theme
proteasome_nodes <- queryMap(cm_obj_adv, query = "BIOCARTA_PROTEASOME_PATHWAY", type = "theme")

cat("Found", nrow(proteasome_nodes), "nodes in the Proteasome theme.")
#> Found 310 nodes in the Proteasome theme.

proteasome_nodes %>%
  select(complexId, primaryFunctionalDomain, proteinCount, betweenness) %>%
  head(10) %>%
  knitr::kable()
```

| complexId   | primaryFunctionalDomain     | proteinCount | betweenness |
|:------------|:----------------------------|-------------:|------------:|
| CpxMap_0414 | BIOCARTA_PROTEASOME_PATHWAY |            8 |   0.1364494 |
| CpxMap_0508 | Unenriched                  |           42 |   0.1112882 |
| CpxMap_0401 | BIOCARTA_CACAM_PATHWAY      |            6 |   0.1003009 |
| CpxMap_0509 | BIOCARTA_PROTEASOME_PATHWAY |          102 |   0.0917965 |
| CpxMap_0505 | BIOCARTA_PROTEASOME_PATHWAY |           40 |   0.0883700 |
| CpxMap_0506 | BIOCARTA_PROTEASOME_PATHWAY |           40 |   0.0766518 |
| CpxMap_0454 | BIOCARTA_PROTEASOME_PATHWAY |            9 |   0.0659081 |
| CpxMap_0182 | BIOCARTA_CACAM_PATHWAY      |            4 |   0.0656916 |
| CpxMap_0501 | BIOCARTA_PROTEASOME_PATHWAY |           20 |   0.0597822 |
| CpxMap_0148 | BIOCARTA_PROTEASOME_PATHWAY |            4 |   0.0572297 |
