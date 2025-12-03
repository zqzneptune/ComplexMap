# Summarize Major Biological Themes in a Complex Map

Identifies and summarizes physical neighborhoods (themes) in the complex
map using community detection.

## Usage

``` r
summarizeThemes(
  complexMapObject,
  method = "louvain",
  add_to_object = TRUE,
  verbose = TRUE
)
```

## Arguments

- complexMapObject:

  A \`ComplexMap\` object.

- method:

  Community detection method ("louvain", "walktrap", "infomap").
  Defaults to "louvain".

- add_to_object:

  Logical. If \`TRUE\` (default), returns the \`ComplexMap\` object with
  \`themeId\` and \`themeLabel\` columns added to the node table. If
  \`FALSE\`, returns a summary \`tibble\` of the themes.

- verbose:

  Logical.

## Value

If \`add_to_object = TRUE\`, a modified \`ComplexMap\` object. If
\`add_to_object = FALSE\`, a \`tibble\` with columns: - \`themeId\`:
Cluster ID. - \`themeLabel\`: The consensus functional label. -
\`themePurity\`: The fraction of nodes matching the primary label
(0-1). - \`nodeCount\`: Number of complexes. - \`edgeCount\`: Number of
internal edges.

## Details

\*\*Systems Biology Rationale:\*\* Since the network layout is driven
primarily by physical composition (alpha=0.75), the communities detected
here represent \*\*Physical Machines\*\* or neighborhoods.

This function characterizes these machines by their functions. A single
physical machine might contain complexes with slightly different
specific labels. To reflect this, the labeling logic is: 1. Identify the
most frequent functional domain in the cluster. 2. Calculate \*\*Theme
Purity\*\* (what 3. If purity is low (\< 50 (e.g., "Function A /
Function B") to indicate multi-functionality.

By default, this function adds the theme assignments directly to the
node table of the \`ComplexMap\` object, making it easy to visualize or
query by theme.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
