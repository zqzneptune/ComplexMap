# 2. Start a Typical Workflow

## Introduction

Protein complexes are the functional machinery of the cell.
High-throughput experimental methods can identify hundreds of putative
protein complexes in a single experiment. The `ComplexMap` package
provides a comprehensive, end-to-end workflow to process, analyze,
annotate, and visualize such datasets.

The version 2.0 workflow is built around a **Physical-First**
philosophy: we prioritize the physical composition (protein identity) of
complexes to define the map’s layout, while using functional enrichment
to provide a “color landscape” that highlights biological diversity
without collapsing distinct sub-complexes into generic blobs.

## Step 1: Loading Data

The `ComplexMap` package includes the `demoComplexes` dataset (622
putative human complexes). We also provide a helper to locate an example
BioCarta gene set (GMT) file.

``` r

# Load example complexes
data("demoComplexes", package = "ComplexMap")

# Access the included example GMT file (BioCarta)
gmtPath <- getExampleGmt()
biocartaGmt <- getGmtFromFile(gmtPath, verbose = FALSE)

# Preview the input
message(sprintf("Loaded %d complexes and %d gene sets.", 
                length(demoComplexes), length(biocartaGmt)))
#> Loaded 622 complexes and 292 gene sets.
```

## Step 2: Running the Main Workflow

The core of the package is the
[`createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
function. This high-level wrapper handles:

1.  **Refinement**: (Optional) Merges highly redundant complexes to
    clean the input.

2.  **Functional Enrichment**: Calculates p-values and **Fold
    Enrichment** for every complex.

3.  **Network Construction**: Links complexes based on shared proteins
    (`alpha=0.75` gives 75% weight to physical overlap).

4.  **Semantic Clustering**: Groups functional terms to assign a
    consistent color palette that maximizes diversity.

5.  **Topology Calculation**: Runs a Fruchterman-Reingold layout and
    calculates centrality metrics (betweenness and degree).

``` r

# Run the entire workflow
# We use a subset of 200 complexes for this vignette demonstration
cmap <- createComplexMap(
  complexList = demoComplexes[1:200],
  gmt = biocartaGmt,
  ifRefineCpx = TRUE,    # Clean up redundancy
  mergeThreshold = 0.90, # Only merge nearly identical complexes
  alpha = 0.75,          # Layout priority: 75% Physical, 25% Functional
  verbose = TRUE
)
#> --- Starting ComplexMap Workflow ---
#> Parameters: similarity='jaccard', alpha=0.75 (Diversity Priority)
#> Refinement: Enabled
#> 
#> Step 1: Refining complex list...
#> 
#> --- Refining Input Complex List (Minimal Merging Strategy) ---
#> Filtered 34 complexes by size. Retaining 166.
#> Identifying merge groups: method='jaccard', threshold >= 0.90...
#> Found 0 redundancy groups. Merging 166 complexes into 166.
#> 
#> --- Refinement Complete ---
#> 
#> Step 2: Running enrichment analysis...
#> Running enrichment for 166 complexes (Cutoff: 0.05)...
#> Annotation complete. Found significant terms for 64 complexes.
#> 
#> Step 3: Building complex network...
#> Building complex network (combined mode, alpha=0.75)...
#> Processing 13695 pairs using 11 cores...
#> Network built: 94 edges retained.
#> 
#> Step 4: Generating node attributes...
#> Generating node attributes (prioritizing functional specificity colors)...
#>     -> Clustering 135 terms using co-occurrence (jaccard)
#> Metric: 'jaccard' with unit: 'log'; comparing: 135 vectors
#>     -> Generating diverse palette: 25 functional domains (Average Linkage).
#> 
#> Step 5: Computing map topology...
#> Computing map topology (layout and centrality)...
#> Topology computation complete.
#> 
#> --- ComplexMap Workflow Complete ---
```

#### The Systems Biology Dashboard

One of the key features of the `ComplexMap` object is its `print`
method, which provides immediate feedback on the “health” of your
landscape.

``` r

# Inspect the resulting object
print(cmap)
#> # ComplexMap Object (Physical-First Layout)
#> # -- Physical Structure: 166 nodes, 94 edges (0.57 edges/node)
#> # -- Functional Landscape:
#> #    * Diversity: 24 distinct functional domains (colors)
#> #    * Coverage:  38.6% of complexes annotated
#> # -- Accessors: `getNodeTable()`, `getEdgeTable()`
#> # -- Explore:   `explore(cmap)` to launch interactive Cytoscape.js viewer.
```

**What to look for:**

- **Diversity**: A high count of distinct functional domains indicates a
  specific, high-resolution map.

- **Coverage**: Shows what percentage of your experimental data could be
  biologically annotated.

## Step 3: Exploring and Querying

You can interact with the map data using standard accessors or the
targeted `queryMap` tool.

#### 3.1 Targeted Querying

Let’s find complexes containing **“SMAD4”** (a central mediator of
TGF-beta signaling).

``` r

# Search for protein members
smad_complexes <- queryMap(cmap, query = "SMAD4", type = "protein")
#> Warning: No matches for 'SMAD4'.

if (nrow(smad_complexes) > 0) {
  smad_complexes %>%
    select(complexId, primaryFunctionalDomain, proteinCount) %>%
    knitr::kable()
} else {
  message("SMAD4 not found in this subset.")
}
#> SMAD4 not found in this subset.
```

#### 3.2 Accessing Tables

The node and edge tables are available as tibbles for custom downstream
analysis (e.g., in `ggplot2`).

``` r

nodes <- getNodeTable(cmap)
edges <- getEdgeTable(cmap)

# Example: Top 3 most central complexes by betweenness
nodes %>%
  select(complexId, primaryFunctionalDomain, betweenness) %>%
  arrange(desc(betweenness)) %>%
  head(3)
#> # A tibble: 3 × 3
#>   complexId   primaryFunctionalDomain     betweenness
#>   <chr>       <chr>                             <dbl>
#> 1 CpxMap_0090 BIOCARTA_SALMONELLA_PATHWAY      0.0434
#> 2 CpxMap_0121 BIOCARTA_HBX_PATHWAY             0.0344
#> 3 CpxMap_0088 BIOCARTA_HBX_PATHWAY             0.0212
```

## Step 4: Interactive Visualization

In `ComplexMap >= 2.0.0`, static R plots have been replaced by a
powerful interactive engine based on **Cytoscape.js**.

To launch the explorer, simply call
[`explore()`](https://zqzneptune.github.io/ComplexMap/reference/explore.md)
on your object. This will open a Shiny application where you can:

1.  **Zoom and Pan** through the physical landscape.

2.  **Filter** by functional domain or protein count.

3.  **Inspect** specific nodes to see their full protein list and
    enrichment scores.

``` r

# Launch the interactive viewer
explore(cmap)
```

## Step 5: Exporting to External Tools

If you prefer to perform advanced network analysis in **Cytoscape
Desktop** or **Gephi**, use `exportNetwork`. It sanitizes complex R data
types into strings and writes standard TSV or GraphML files.

``` r

# Export for Cytoscape Desktop
exportNetwork(cmap, filePrefix = "my_complex_landscape", format = "cytoscape")
```

## Conclusion

The `ComplexMap` workflow transforms a raw list of protein
identifications into a structured, functional landscape. By prioritizing
physical composition and enforcing functional diversity, it ensures that
your map reflects the true biological complexity of the proteome.
