# Export a ComplexMap Object to Cytoscape.js JSON

Converts a \`ComplexMap\` object into a Cytoscape.js-compatible JSON
string, with scaled layout positions. This JSON is consumed directly by
the Shiny/Cytoscape.js interactive explorer.

## Usage

``` r
toCytoscapeJSON(
  cmap,
  prune = TRUE,
  prune_method = "top_k",
  k = 5,
  weight_quantile = 0.75,
  coord_range = c(-500, 500),
  verbose = TRUE
)
```

## Arguments

- cmap:

  A \`ComplexMap\` object (with layout coordinates in the node table).

- prune:

  Logical. Whether to prune the network before export. Default \`TRUE\`.

- prune_method:

  Pruning method passed to \`pruneNetwork()\`. Default \`"top_k"\`.

- k:

  Integer. Top-k edges per node if \`prune_method = "top_k"\`. Default
  \`5\`.

- weight_quantile:

  Numeric. Quantile threshold if \`prune_method = "quantile"\`. Default
  \`0.75\`.

- coord_range:

  Numeric length-2 vector for coordinate scaling. Default \`c(-500,
  500)\`.

- verbose:

  Logical.

## Value

A JSON string suitable for \`window.complexmapElements\` in JavaScript.
