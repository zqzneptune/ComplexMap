# Visualize a Complex Map with Direct Node Labels

Creates a static visualization of the complex network using \`ggraph\`,
where functional domain labels are placed directly on the plot near the
nodes.

## Usage

``` r
visualizeMapDirectLabels(
  layoutDf,
  edgesDf,
  title = "ComplexMap Functional Landscape",
  subtitle = "Nodes are protein complexes, colored by function",
  bgColor = "black",
  edgeColor = "white",
  nodeSizeRange = c(2, 10),
  labelFillColor = ggplot2::alpha("white", 0.7),
  fontFamily = "sans",
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

- labelFillColor:

  The background fill color for the labels.

- fontFamily:

  The base font family for all plot text.

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A \`ggplot\` object representing the network visualization.

## Details

This function is designed for clarity when domains are distinct. For
large, dense domains (more than 2 complexes), a single label is placed
at the domain's centroid. For smaller domains, each complex is labeled
individually. It uses the pre-computed layout from
\`computeMapTopology\`.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
nodes <- tibble::tibble(
  complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
  primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Metabolism"),
  sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#0000FF")
)
edges <- tibble::tibble(
  source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
)

# --- Generate Plot ---
if (requireNamespace("ggrepel", quietly = TRUE)) {
  visualizeMapDirectLabels(nodes, edges)
}
#> Visualizing ComplexMap with direct labels...

```
