# Visualize a Complex Map with Direct Node Labels

Creates a static visualization where functional labels are placed
directly on the plot.

## Usage

``` r
visualizeMapDirectLabels(
  layoutDf,
  edgesDf,
  title = "ComplexMap Functional Landscape",
  subtitle = "Nodes colored by Specificity-Weighted Function",
  bgColor = "black",
  edgeColor = "white",
  nodeSizeRange = c(2, 10),
  labelFillColor = ggplot2::alpha("white", 0.7),
  centroid_threshold = 2,
  color.by = NULL,
  color.palette = "viridis",
  color.legend.title = NULL,
  verbose = TRUE
)
```

## Arguments

- layoutDf:

  A data frame containing node attributes and layout coordinates.

- edgesDf:

  A data frame containing the network edges.

- title:

  Plot title.

- subtitle:

  Plot subtitle.

- bgColor:

  Background color.

- edgeColor:

  Edge color.

- nodeSizeRange:

  Min and max node size.

- labelFillColor:

  Background fill color for labels.

- centroid_threshold:

  Domains with \> N complexes get a single centroid label.

- color.by:

  Optional numeric column for continuous coloring.

- color.palette:

  Palette name ("viridis", "plasma") or vector of colors.

- color.legend.title:

  Title for color legend.

- verbose:

  Logical.

## Value

A \`ggplot\` object.
