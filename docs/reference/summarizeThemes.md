# Summarize Major Biological Themes in a Complex Map

This function analyzes the network structure within a \`ComplexMap\`
object to identify and summarize major biological themes. It uses
community detection algorithms to find densely connected clusters of
nodes (themes).

## Usage

``` r
summarizeThemes(complexMapObject, method = "louvain", verbose = TRUE)
```

## Arguments

- complexMapObject:

  A \`ComplexMap\` object returned by \`createComplexMap()\`.

- method:

  A character string specifying the community detection algorithm to
  use. Must be a valid \`igraph\` clustering function (e.g., "louvain",
  "walktrap", "infomap"). Defaults to "louvain".

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A \`tibble\` where each row represents a summarized theme. The tibble
contains the following columns: - \`themeId\`: A unique integer
identifier for the theme. - \`themeLabel\`: A descriptive label for the
theme. - \`nodeCount\`: The number of nodes (complexes) in the theme. -
\`edgeCount\`: The number of internal edges within the theme.

## Details

The function performs the following steps: 1. It constructs an
\`igraph\` graph object from the node and edge tables of the
\`complexMapObject\`. 2. It applies a community detection algorithm
(e.g., Louvain, as default) to partition the network into clusters or
"themes". 3. For each identified theme, it generates a descriptive
\`themeLabel\` by finding the most frequently occurring
\`primaryFunctionalDomain\` among the member nodes (excluding
"Unenriched"). 4. It calculates summary statistics for each theme,
including the number of nodes and edges it contains.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Assume 'cm_obj' is a valid ComplexMap object created by createComplexMap()
# if (requireNamespace("igraph", quietly = TRUE)) {
#   themeSummary <- summarizeThemes(cm_obj)
#   print(themeSummary)
# }
```
