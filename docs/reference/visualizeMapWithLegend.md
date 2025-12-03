# Visualize a Complex Map with a Color Legend

Creates a static visualization with a discrete legend for functional
domains.

## Usage

``` r
visualizeMapWithLegend(
  layoutDf,
  edgesDf,
  title = "ComplexMap Functional Landscape",
  subtitle = "Nodes colored by Specificity-Weighted Function",
  bgColor = "black",
  edgeColor = "white",
  nodeSizeRange = c(2, 10),
  unenrichedColor = "#CCCCCC",
  color.by = NULL,
  color.palette = "plasma",
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

- unenrichedColor:

  Color for unenriched nodes.

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
