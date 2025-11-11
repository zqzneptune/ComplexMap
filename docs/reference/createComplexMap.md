# Create a Complete Complex Map Object

A high-level wrapper function that executes the entire \`ComplexMap\`
workflow. It takes a list of protein complexes and a gene set matrix
(GMT) and performs refinement, enrichment, network construction, and
topology calculation.

## Usage

``` r
createComplexMap(complexList, gmt, verbose = TRUE, ...)
```

## Arguments

- complexList:

  A named list where each element is a character vector of protein
  identifiers representing a complex.

- gmt:

  A named list where each element is a character vector of genes,
  representing a functional gene set (e.g., from a GMT file).

- verbose:

  A logical value indicating whether to print progress messages for the
  entire workflow. Defaults to \`TRUE\`.

- ...:

  Additional arguments to be passed down to the core functions. Common
  arguments include: - \`minSize\`, \`maxSize\`, \`mergeThreshold\` (for
  \`refineComplexList\`) - \`pAdjustMethod\`, \`pValueCutoff\` (for
  \`runComplexEnrichment\`) - \`mode\`, \`similarityMethod\`, \`alpha\`
  (for \`buildComplexNetwork\`)

## Value

A validated \`ComplexMap\` S3 object containing the final node and edge
tables.

## Details

This function serves as the primary entry point for most analyses. It
internally calls the core workflow in the following order: 1.
\`refineComplexList()\` 2. \`runComplexEnrichment()\` 3.
\`buildComplexNetwork()\` 4. \`generateNodeAttributes()\` 5.
\`computeMapTopology()\`

Arguments for the underlying functions can be passed directly to this
wrapper via the \`...\` parameter.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Assume 'demoComplexes' and a 'gmt' object are loaded
# gmtPath <- getExampleGmt()
# gmt <- getGmtFromFile(gmtPath, verbose = FALSE)

# Run the full workflow with custom parameters
# complexMapObject <- createComplexMap(
#   demoComplexes,
#   gmt,
#   verbose = TRUE,
#   minSize = 5,
#   mergeThreshold = 0.8,
#   pValueCutoff = 0.01,
#   mode = "combined"
# )
# print(complexMapObject)
```
