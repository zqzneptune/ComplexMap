# Prune Network Edges Before Cytoscape Export

Reduces the number of edges in a \`ComplexMap\` object to prevent
browser overload during interactive exploration.

## Usage

``` r
pruneNetwork(
  cmap,
  method = "top_k",
  k = 5,
  weight_quantile = 0.75,
  verbose = TRUE
)
```

## Arguments

- cmap:

  A \`ComplexMap\` object.

- method:

  Pruning method: \`"top_k"\` (keep top-k edges per node) or
  \`"quantile"\` (keep edges above a weight quantile).

- k:

  Integer. Number of top edges to retain per node when \`method =
  "top_k"\`.

- weight_quantile:

  Numeric (0-1). Quantile threshold when \`method = "quantile"\`.
  Defaults to \`0.75\`.

- verbose:

  Logical.

## Value

A modified \`ComplexMap\` object with pruned edges.
