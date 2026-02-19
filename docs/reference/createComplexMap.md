# Create a Complete Complex Map Object

A high-level wrapper function that executes the entire \`ComplexMap\`
workflow. It takes a list of protein complexes and a gene set matrix
(GMT) and performs refinement, enrichment, network construction, and
topology calculation.

## Usage

``` r
createComplexMap(
  complexList,
  gmt,
  similarityMethod = "jaccard",
  alpha = 0.75,
  ifRefineCpx = FALSE,
  layout_seed = 123,
  verbose = TRUE,
  ...
)
```

## Arguments

- complexList:

  A named list where each element is a character vector of protein
  identifiers representing a complex.

- gmt:

  A named list where each element is a character vector of genes.

- similarityMethod:

  The metric used for comparing complexes and functional terms. Defaults
  to \*\*"jaccard"\*\* to penalize size differences and maintain
  diversity. Avoid "overlap" if you wish to see specific pathways
  distinct from generic ones.

- alpha:

  Numeric (0-1). The weight given to physical protein composition versus
  functional similarity in the network layout. Defaults to \*\*0.75\*\*
  (75% physical, 25% functional).

- ifRefineCpx:

  Logical. If \*\*TRUE\*\*, executes Step 1 to refine the input complex
  list (merge redundant complexes). If \*\*FALSE\*\*, Step 1 is skipped,
  and the raw \`complexList\` is used for downstream analysis. Defaults
  to \*\*FALSE\*\*.

- layout_seed:

  Integer seed for reproducible layout. Defaults to 123.

- verbose:

  A logical value indicating whether to print progress messages.

- ...:

  Additional arguments passed to core functions: - \`minSize\`,
  \`maxSize\`, \`mergeThreshold\` (for \`refineComplexList\` if
  \`ifRefineCpx = TRUE\`) - \`pAdjustMethod\`, \`pValueCutoff\` (for
  \`runComplexEnrichment\`) - \`geneSetDb\` (for
  \`generateNodeAttributes\` semantic clustering)

## Value

A validated \`ComplexMap\` S3 object.

## Details

This function is tuned to generate a \*\*functionally diverse
landscape\*\*. It enforces the following logic: 1. \*\*Refinement:\*\*
(Optional) Uses Jaccard similarity to merge only highly redundant
complexes, preserving biological variants (subsets/supersets).
Controlled by \`ifRefineCpx\`. 2. \*\*Enrichment:\*\* Calculates 'Fold
Enrichment' to prioritize specific biological functions over generic
ones. 3. \*\*Network:\*\* Builds the layout primarily based on
\*\*Physical Composition\*\* (alpha = 0.75), using functional
annotations only to group related clusters, not to collapse them.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
