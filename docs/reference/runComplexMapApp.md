# Run the ComplexMap Interactive Explorer

Launches a Shiny application that renders a \`ComplexMap\` object as an
interactive Cytoscape.js network. Nodes are colored by functional
domain, sized by protein count, and fully filterable via sidebar
controls.

## Usage

``` r
runComplexMapApp(
  cmap,
  prune = TRUE,
  prune_method = "top_k",
  k = 5,
  weight_quantile = 0.75,
  port = NULL,
  launch.browser = TRUE,
  verbose = TRUE
)
```

## Arguments

- cmap:

  A \`ComplexMap\` object.

- prune:

  Logical. Whether to prune edges before launching. Default \`TRUE\`.

- prune_method:

  Pruning method: \`"top_k"\` or \`"quantile"\`. Default \`"top_k"\`.

- k:

  Integer. Top-k edges per node when \`prune_method = "top_k"\`. Default
  \`5\`.

- weight_quantile:

  Numeric (0-1). Quantile threshold for \`"quantile"\` method. Default
  \`0.75\`.

- port:

  Optional integer port for the Shiny server.

- launch.browser:

  Logical. Open browser automatically. Default \`TRUE\`.

- verbose:

  Logical.

## Value

Launches the Shiny app (invisibly returns \`NULL\`).
