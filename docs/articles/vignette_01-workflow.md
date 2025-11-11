# Analyzing a Human Protein Complex Map

## Introduction

Protein complexes are the functional machinery of the cell.
High-throughput experimental methods, such as co-fractionation followed
by mass spectrometry (CF-MS), can identify hundreds of putative protein
complexes in a single experiment. The `ComplexMap` package provides a
comprehensive workflow to process, analyze, annotate, and visualize such
a dataset.

This vignette demonstrates a typical workflow using a dataset of human
soluble protein complexes identified via co-fractionation, originally
published in 2012 [A census of human soluble protein
complexes](https://doi.org/10.1016/j.cell.2012.08.011). We will perform
quality control, refine the list, perform functional enrichment, and
build a network-based “map” of the functional landscape of these
complexes.

First, let’s load the `ComplexMap` package, along with `dplyr` for data
manipulation.

``` r
library(ComplexMap)
library(dplyr)
```

## Step 1: Loading Data

The `ComplexMap` package includes two key example datasets:
`demoComplexes` and `referenceComplexes`.

- `demoComplexes`: A dataset containing a list of 622 putative protein
  complexes identified in our study using an integrative global
  proteomic profiling approach.

- `referenceComplexes`: A reference dataset containing 324 merged CORUM
  protein complexes used in training protein-protein interaction scoring
  and clustering optimization procedures in the study. These complexes
  were curated from the CORUM database and merged to reduce redundancy
  for benchmarking purposes.

We load these datasets directly.

``` r
# Load the example datasets shipped with the package
data("demoComplexes")
data("referenceComplexes")
```

## Step 2: Quality Control

Before any analysis, it’s crucial to assess the quality of the input
complex list. The
[`qcComplexList()`](https://zqzneptune.github.io/ComplexMap/reference/qcComplexList.md)
function provides a summary of basic statistics, complex sizes, and
pairwise redundancy.

``` r
qcComplexList(demoComplexes)
```

    ## 
    ## --- Running Quality Control on Complex List ---

    ## 
    ## [1] Basic Statistics:

    ##     - Total number of complexes: 622

    ##     - Total number of unique proteins: 2649

    ## 
    ## [2] Complex Size Distribution:

    ##        Min. 1st Qu.  Median    Mean 3rd Qu.    Max.

    ##       1.000   3.000   4.000   5.183   5.000 102.000

    ## Warning in qcComplexList(demoComplexes): 112 complexes have fewer than 3
    ## members.

    ## 
    ## [3] Redundancy Analysis:

    ##     - Distribution of Jaccard similarity scores:

    ##          Min.   1st Qu.    Median      Mean   3rd Qu.      Max.

    ##     0.0000000 0.0000000 0.0000000 0.0003226 0.0000000 0.5000000

    ##     - No highly redundant complex pairs detected.

    ## 
    ## --- QC Complete ---

The QC report tells us we have 622 complexes composed of 2649 unique
proteins. It also warns if there some complexes are highly redundant
(Jaccard similarity \>= 0.8), suggesting that our next step, refinement,
is necessary.

## Step 3: Refining the Complex List

The initial list of complexes may contain very small complexes or highly
overlapping ones that represent slight variations of the same biological
entity. The
[`refineComplexList()`](https://zqzneptune.github.io/ComplexMap/reference/refineComplexList.md)
function addresses this by filtering by size and merging redundant
complexes.

Here, we will merge any complexes that have a Jaccard similarity of 0.75
or higher.

``` r
# We use a slightly lower mergeThreshold to be more aggressive for this demo
refinedComplexes <- refineComplexList(demoComplexes, mergeThreshold = 0.75)
```

    ## 
    ## --- Refining Input Complex List ---

    ## Filtered 112 complexes by size. Retaining 510.

    ## Identifying merge groups with Jaccard >= 0.75...

    ## Found 0 merge groups. Merging 0 complexes into 510.

    ## Merging complete. Final list has 510 complexes.

    ## 
    ## --- Refinement Complete ---

``` r
# Let's see how many complexes we have now
length(refinedComplexes)
```

    ## [1] 510

After refinement, the list has been consolidated to 510 complexes, which
provides a cleaner basis for downstream analysis.

## Step 4: Functional Enrichment

To understand the biological role of each complex, we can perform a
functional enrichment analysis. We first need a gene set database, such
as the BioCarta pathways. The package includes a helper function to
access an example GMT file.

``` r
# Get the path to the example GMT file
gmtPath <- getExampleGmt()
biocartaGmt <- getGmtFromFile(gmtPath)
```

    ## Fetching gene sets from local file: /tmp/RtmpXXrxdv/temp_libpath31341c885c94/ComplexMap/extdata/c2.cp.biocarta.v2025.1.Hs.symbols.gmt

``` r
# Run enrichment analysis on the refined complex list
enrichments <- runComplexEnrichment(refinedComplexes, biocartaGmt)
```

    ## Running enrichment for 510 complexes...

    ## Annotation complete. Found terms for 239 complexes.

``` r
# View the enrichment results for the first complex with significant terms
if (length(enrichments) > 0) {
  head(enrichments[])
}
```

    ## $CpxMap_0006
    ##                            ID                 Description   p.adjust Count
    ## 1 BIOCARTA_EICOSANOID_PATHWAY BIOCARTA_EICOSANOID_PATHWAY 0.01457919     1
    ##   Ratio     Fold
    ## 1     1 68.59091
    ## 
    ## $CpxMap_0008
    ##                         ID              Description   p.adjust Count Ratio
    ## 1     BIOCARTA_EPO_PATHWAY     BIOCARTA_EPO_PATHWAY 0.02208968     1     1
    ## 2      BIOCARTA_GH_PATHWAY      BIOCARTA_GH_PATHWAY 0.02236581     1     1
    ## 3   BIOCARTA_IL2RB_PATHWAY   BIOCARTA_IL2RB_PATHWAY 0.02451955     1     1
    ## 4     BIOCARTA_IL3_PATHWAY     BIOCARTA_IL3_PATHWAY 0.02208968     1     1
    ## 5 BIOCARTA_NKCELLS_PATHWAY BIOCARTA_NKCELLS_PATHWAY 0.02208968     1     1
    ##        Fold
    ## 1  79.42105
    ## 2  55.88889
    ## 3  40.78378
    ## 4 100.60000
    ## 5  75.45000
    ## 
    ## $CpxMap_0011
    ##                     ID          Description    p.adjust Count Ratio     Fold
    ## 1 BIOCARTA_RAB_PATHWAY BIOCARTA_RAB_PATHWAY 0.008614977     1     1 116.0769
    ## 
    ## $CpxMap_0012
    ##                    ID         Description   p.adjust Count Ratio   Fold
    ## 1 BIOCARTA_G2_PATHWAY BIOCARTA_G2_PATHWAY 0.01590457     1     1 62.875
    ## 
    ## $CpxMap_0017
    ##                              ID                   Description   p.adjust Count
    ## 1      BIOCARTA_CASPASE_PATHWAY      BIOCARTA_CASPASE_PATHWAY 0.04430150     1
    ## 2        BIOCARTA_DEATH_PATHWAY        BIOCARTA_DEATH_PATHWAY 0.04430150     1
    ## 3  BIOCARTA_DNAFRAGMENT_PATHWAY  BIOCARTA_DNAFRAGMENT_PATHWAY 0.03964278     1
    ## 4          BIOCARTA_FAS_PATHWAY          BIOCARTA_FAS_PATHWAY 0.04430150     1
    ## 5 BIOCARTA_MITOCHONDRIA_PATHWAY BIOCARTA_MITOCHONDRIA_PATHWAY 0.04430150     1
    ## 6       BIOCARTA_PARKIN_PATHWAY       BIOCARTA_PARKIN_PATHWAY 0.03964278     1
    ## 7          BIOCARTA_SET_PATHWAY          BIOCARTA_SET_PATHWAY 0.03964278     1
    ## 8        BIOCARTA_TNFR1_PATHWAY        BIOCARTA_TNFR1_PATHWAY 0.04430150     1
    ##   Ratio     Fold
    ## 1   0.5 34.29545
    ## 2   0.5 26.01724
    ## 3   0.5 75.45000
    ## 4   0.5 25.15000
    ## 5   0.5 39.71053
    ## 6   0.5 94.31250
    ## 7   0.5 75.45000
    ## 8   0.5 26.01724
    ## 
    ## $CpxMap_0020
    ##                     ID          Description   p.adjust Count Ratio     Fold
    ## 1 BIOCARTA_TEL_PATHWAY BIOCARTA_TEL_PATHWAY 0.01126574     1     1 88.76471

The analysis found significant terms for 239 of our refined complexes.
The output is a list where each element is a data frame of significant
functional terms for the corresponding complex.

## Step 5: Building the Complex Network

To visualize the relationships between complexes, we build a similarity
network. An edge is drawn between two complexes if they are similar,
either by sharing proteins (compositional similarity) or by sharing
functions (functional similarity).

The
[`buildComplexNetwork()`](https://zqzneptune.github.io/ComplexMap/reference/buildComplexNetwork.md)
function calculates these similarities and creates an edge list. We’ll
use a “combined” mode, which creates a weighted average of both
similarity types.

``` r
networkEdges <- buildComplexNetwork(
  complexes = refinedComplexes,
  enrichments = enrichments,
  mode = "combined",
  similarityMethod = "jaccard"
)
```

    ## Building complex network using 'jaccard' similarity...

    ## Using 11 cores for parallel processing.

    ## Processing 129795 complex pairs...

    ## Split into 130 chunks of up to 1000 pairs each.

    ## Combining results from chunks...

    ## Calculating final weights and filtering...

    ## Network construction complete: 1390 edges retained.

``` r
glimpse(networkEdges)
```

    ## Rows: 1,390
    ## Columns: 8
    ## $ source_complex_id <chr> "CpxMap_0005", "CpxMap_0006", "CpxMap_0008", "CpxMap…
    ## $ target_complex_id <chr> "CpxMap_0453", "CpxMap_0356", "CpxMap_0090", "CpxMap…
    ## $ compSim           <dbl> 0.09090909, 0.00000000, 0.00000000, 0.00000000, 0.00…
    ## $ funcSim           <dbl> 0.00000000, 0.05555556, 0.04166667, 0.20000000, 0.36…
    ## $ sharedProt        <int> 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0…
    ## $ sharedFunc        <int> 0, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    ## $ weight            <dbl> 0.04545455, 0.02777778, 0.02083333, 0.10000000, 0.18…
    ## $ similarity_mode   <chr> "combined", "combined", "combined", "combined", "com…

The resulting `tibble` contains 1390 edges, each with calculated
similarity scores and a final `weight` that will be used for
visualization.

## Step 6: Generating Node Attributes for Visualization

Before plotting, we need to generate attributes for each node (complex)
in the network, such as size, color, and a primary functional label. The
[`generateNodeAttributes()`](https://zqzneptune.github.io/ComplexMap/reference/generateNodeAttributes.md)
function is designed for this. It clusters enriched terms into
“functional domains” and calculates a unique “blended” color for each
complex based on its functional profile.

``` r
nodeAttributes <- generateNodeAttributes(
  complexes = refinedComplexes,
  enrichments = enrichments
)
```

    ## Generating core node attributes (function and color)...

    ##     -> Clustering terms using 'jaccard' similarity.

    ## Metric: 'jaccard' with unit: 'log'; comparing: 231 vectors

``` r
glimpse(nodeAttributes)
```

    ## Rows: 510
    ## Columns: 7
    ## $ complexId               <chr> "CpxMap_0001", "CpxMap_0002", "CpxMap_0003", "…
    ## $ proteinCount            <int> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3…
    ## $ proteins                <chr> "WDR82,ULK3,UCHL1", "PPP1R8,PPP1CC,NF2", "IFT8…
    ## $ primaryFunctionalDomain <chr> "Unenriched", "Unenriched", "Unenriched", "Une…
    ## $ topEnrichedFunctions    <chr> NA, NA, NA, NA, NA, "BIOCARTA_EICOSANOID_PATHW…
    ## $ colorHex                <chr> "#CCCCCC", "#CCCCCC", "#CCCCCC", "#CCCCCC", "#…
    ## $ sizeMapping             <dbl> 1.584963, 1.584963, 1.584963, 1.584963, 1.5849…

## Step 7: Computing Network Topology

With the nodes and edges defined, the final step before plotting is to
calculate the network layout. The
[`computeMapTopology()`](https://zqzneptune.github.io/ComplexMap/reference/computeMapTopology.md)
function uses a force-directed algorithm to determine the (x, y)
coordinates for each node and calculates centrality metrics like degree
and betweenness.

``` r
mapLayout <- computeMapTopology(nodeAttributes, networkEdges)
```

    ## Computing map topology (layout and centrality)...

    ## Topology computation complete.

``` r
glimpse(mapLayout)
```

    ## Rows: 510
    ## Columns: 11
    ## $ complexId               <chr> "CpxMap_0359", "CpxMap_0414", "CpxMap_0508", "…
    ## $ proteinCount            <int> 5, 8, 42, 6, 102, 4, 40, 40, 4, 9, 4, 20, 15, …
    ## $ proteins                <chr> "PRKACB,PRKACA,PRKAR2A,CIRBP,CAPRIN1", "WDR3,R…
    ## $ primaryFunctionalDomain <chr> "BIOCARTA_EICOSANOID_PATHWAY", "BIOCARTA_EICOS…
    ## $ topEnrichedFunctions    <chr> "BIOCARTA_AGPCR_PATHWAY, BIOCARTA_AKAP13_PATHW…
    ## $ colorHex                <chr> "#E41A1C", "#E41A1C", "#CCCCCC", "#E41A1C", "#…
    ## $ sizeMapping             <dbl> 2.321928, 3.000000, 5.392317, 2.584963, 6.6724…
    ## $ x                       <dbl> 7.9898939, 6.7094059, -3.2063886, 10.0283948, …
    ## $ y                       <dbl> 6.1280771, 5.0562320, 2.5393381, 4.6364274, -0…
    ## $ betweenness             <dbl> 0.15447148, 0.13644942, 0.11128815, 0.10030088…
    ## $ degree                  <dbl> 33, 45, 20, 28, 22, 41, 19, 7, 26, 13, 22, 27,…

The `mapLayout` data frame is now the “master” table, containing all the
information needed for plotting.

## Step 8: Visualization

`ComplexMap` provides three ways to visualize the final network.

#### 8.1 Static Plot with Direct Labels

This plot is excellent for publications, where labels are placed
directly onto the plot.

``` r
# ggrepel is required for this plot
if (requireNamespace("ggrepel", quietly = TRUE)) {
  visualizeMapDirectLabels(mapLayout, networkEdges)
}
```

    ## Visualizing ComplexMap with direct labels...

![A network plot of protein complexes with functional domain labels
placed directly on the
map.](vignette_01-workflow_files/figure-html/vis-direct-labels-1.png)

#### 8.2 Static Plot with a Legend

This version is useful when direct labels are too cluttered. It uses a
discrete color legend to represent the functional domains.

``` r
visualizeMapWithLegend(mapLayout, networkEdges)
```

    ## Visualizing ComplexMap with a color legend...

![A network plot of protein complexes where node color corresponds to a
functional domain listed in a
legend.](vignette_01-workflow_files/figure-html/vis-legend-1.png)

#### 8.3 Interactive Plot

For deep exploration, an interactive HTML widget is ideal. You can zoom,
pan, and hover over nodes to see detailed tooltips.

``` r
# visNetwork is required for this plot
if (requireNamespace("visNetwork", quietly = TRUE)) {
  visualizeMapInteractive(mapLayout, networkEdges)
}
```

    ## Generating interactive visNetwork plot...

## Conclusion

This vignette has demonstrated the full workflow of the `ComplexMap`
package, from raw complex list to insightful, publication-quality
visualizations. By integrating quality control, refinement, enrichment,
and network analysis, `ComplexMap` provides a powerful and user-friendly
platform for exploring the landscape of protein complexes.
