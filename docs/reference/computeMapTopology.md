# Compute Network Topology and Layout Coordinates

Calculates the 2D layout coordinates and key centrality metrics for a
complex network. This function serves as the final step in preparing
network data for visualization.

## Usage

``` r
computeMapTopology(nodeAttributes, network, verbose = TRUE)
```

## Arguments

- nodeAttributes:

  A \`tibble\` or \`data.frame\` containing attributes for each node
  (complex). Must contain a column with node identifiers that matches
  the source/target columns in the \`network\` data.

- network:

  A \`tibble\` or \`data.frame\` representing the network edges. Must
  contain columns for source, target, and edge \`weight\`.

- verbose:

  A logical value indicating whether to print progress messages.
  Defaults to \`TRUE\`.

## Value

A \`tibble\` containing all original columns from \`nodeAttributes\`
plus four new columns: \`x\`, \`y\` (layout coordinates),
\`betweenness\`, and \`degree\`. The table is arranged in descending
order of betweenness and degree.

## Details

This function takes a node attribute table and an edge list (network)
and performs the following steps: 1. Constructs an \`igraph\` graph
object from the provided data. 2. Computes a force-directed layout using
the Fruchterman-Reingold algorithm via \`ggraph::create_layout\`. Edge
weights are used to influence the layout, pulling strongly connected
nodes closer together. 3. Calculates node centrality metrics: -
\*\*Betweenness Centrality:\*\* Measures how often a node lies on the
shortest path between other nodes (normalized). - \*\*Degree
Centrality:\*\* The number of edges connected to a node. 4. Merges the
layout coordinates and centrality scores back into the original node
attribute table.

## See also

\`generateNodeAttributes()\`, \`buildComplexNetwork()\`

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
# 1. Node attributes
nodes <- tibble::tibble(
  complexId = c("Cpx1", "Cpx2", "Cpx3"),
  proteinCount = c(10, 8, 12)
)

# 2. Network edges
net <- tibble::tibble(
  source = c("Cpx1", "Cpx2"),
  target = c("Cpx2", "Cpx3"),
  weight = c(0.8, 0.6)
)

# --- Compute Topology ---
masterLayout <- computeMapTopology(nodes, net)
#> Computing map topology (layout and centrality)...
#> Topology computation complete.
print(masterLayout)
#> # A tibble: 3 × 6
#>   complexId proteinCount     x        y betweenness degree
#>   <chr>            <dbl> <dbl>    <dbl>       <dbl>  <dbl>
#> 1 Cpx2                 8 -3.25  0.00972           1      2
#> 2 Cpx1                10 -2.42  0.925             0      1
#> 3 Cpx3                12 -4.16 -1.00              0      1
```
