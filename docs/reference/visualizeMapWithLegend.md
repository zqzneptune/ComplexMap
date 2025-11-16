# Visualize a Complex Map with a Color Legend

Creates a static visualization of the complex network using \`ggraph\`,
with a legend representing node colors.

## Usage

``` r
visualizeMapWithLegend(
  layoutDf,
  edgesDf,
  title = "ComplexMap Functional Landscape",
  subtitle = "Nodes are protein complexes, colored by function",
  bgColor = "black",
  edgeColor = "white",
  nodeSizeRange = c(2, 10),
  unenrichedColor = "#CCCCCC",
  fontFamily = "sans",
  size.legend.title = "Complex Size (log2)",
  color.by = NULL,
  color.palette = "viridis",
  color.legend.title = NULL,
  verbose = TRUE
)
```

## Arguments

- layoutDf:

  A data frame containing node attributes and layout coordinates,
  typically from \`computeMapTopology\`.

- edgesDf:

  A data frame containing the network edges.

- title:

  The main title for the plot.

- subtitle:

  The subtitle for the plot.

- bgColor:

  The background color of the plot.

- edgeColor:

  The color of the network edges.

- nodeSizeRange:

  A numeric vector of length 2 specifying the min and max node size.

- unenrichedColor:

  The color for nodes in the "Unenriched" category in categorical mode.

- fontFamily:

  The base font family for all plot text.

- size.legend.title:

  The title for the node size legend.

- color.by:

  A character string specifying the name of a numeric column in
  \`layoutDf\` to use for continuous node coloring. If \`NULL\`
  (default), the categorical \`colorHex\` column is used.

- color.palette:

  A character string or vector of colors for the continuous gradient
  (e.g., "viridis", "plasma", or \`c("blue", "white", "red")\`). Only
  used when \`color.by\` is specified.

- color.legend.title:

  A character string for the title of the continuous color legend.
  Defaults to the value of \`color.by\`.

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A \`ggplot\` object representing the network visualization.

## Details

This function supports two coloring modes:

1\. \*\*Categorical (default):\*\* When \`color.by = NULL\`, nodes are
colored by their \`primaryFunctionalDomain\`. A discrete color legend is
shown.

2\. \*\*Continuous:\*\* When \`color.by\` is set to the name of a
numeric column in \`layoutDf\`, nodes are colored along a continuous
gradient. A color bar legend is shown.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
nodes <- tibble::tibble(
  complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
  primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Unenriched"),
  sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#CCCCCC"),
  purity = c(0.95, 0.87, 0.91)
)
edges <- tibble::tibble(
  source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
)

# --- Usage 1: Default categorical coloring ---
visualizeMapWithLegend(nodes, edges)
#> Visualizing ComplexMap with a color legend...


# --- Usage 2: Continuous coloring ---
visualizeMapWithLegend(nodes, edges, color.by = "purity")
#> Visualizing ComplexMap with a color legend...

```
