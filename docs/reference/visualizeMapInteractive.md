# Visualize a Complex Map Interactively

Creates a dynamic \`visNetwork\` widget.

## Usage

``` r
visualizeMapInteractive(
  layoutDf,
  edgesDf,
  width = "100%",
  height = "90vh",
  title = "ComplexMap Functional Landscape",
  physicsEnabled = FALSE,
  color.by = NULL,
  color.palette = "viridis",
  verbose = TRUE
)
```

## Arguments

- layoutDf:

  A data frame containing node attributes and layout coordinates.

- edgesDf:

  A data frame containing the network edges.

- width:

  Widget width.

- height:

  Widget height.

- title:

  Plot title.

- physicsEnabled:

  Enable physics simulation.

- color.by:

  Optional numeric column for continuous coloring.

- color.palette:

  Palette name ("viridis", "plasma") or vector of colors.

- verbose:

  Logical.

## Value

A \`visNetwork\` widget.
